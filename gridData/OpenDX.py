# gridData --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.

r"""
:mod:`OpenDX` --- routines to read and write simple OpenDX files
================================================================

The OpenDX format for multi-dimensional grid data. OpenDX is a free
visualization software, see http://www.opendx.org.

.. Note:: This module only implements a primitive subset, sufficient
          to represent n-dimensional regular grids.

The OpenDX scalar file format is specified in Appendix `B.2 Data
Explorer Native Files`_ [#OpenDXformat]_.

If you want to build a dx object from your data you can either use the
convenient :class:`~gridData.core.Grid` class from the top level
module (:class:`gridData.Grid`) or see the lower-level methods
described below.


Building a dx object from a numpy array ``A``
---------------------------------------------

If you have a numpy array ``A`` that represents a density in cartesian
space then you can construct a dx object (named a *field* in OpenDX
parlance) if you provide some additional information that fixes the
coordinate system in space and defines the units along the axes.

The following data are required:

grid
    numpy nD array (typically a nD histogram)
grid.shape
    the shape of the array
origin
    the cartesian coordinates of the center of the (0,0,..,0) grid cell
delta
    :math:`n \times n` array with the length of a grid cell along
    each axis; for regular rectangular grids the off-diagonal
    elements are 0 and the diagonal ones correspond to the
    'bin width' of the histogram, eg ``delta[0,0] = 1.0`` (Angstrom)

For example, to build a :class:`field`::

  dx = OpenDX.field('density')
  dx.add('positions', OpenDX.gridpositions(1, grid.shape, origin, delta))
  dx.add('connections', OpenDX.gridconnections(2, grid.shape))
  dx.add('data', OpenDX.array(3, grid))

or all with the constructor::

  dx = OpenDX.field('density', components=dict(
            positions=OpenDX.gridpositions(1,grid.shape, d.origin, d.delta),
            connections=OpenDX.gridconnections(2, grid.shape),
            data=OpenDX.array(3, grid)))


Building a dx object from a dx file
-----------------------------------

One can also read data from an existing dx file::

 dx = OpenDX.field(0)
 dx.read('file.dx')

The dx :class:`field` object has a method
:meth:`~OpenDX.field.histogramdd` that produces output identical to the
:func:`numpy.histogramdd` function. In this way, one can store nD
histograms in a portable and universal manner::

  histogram, edges = dx.histogramdd()

.. rubric:; Footnotes

.. [#OpenDXformat] The original link to the OpenDX file format specs
   http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF is dead so I am linking
   to an archived copy at the Internet Archive , `B.2 Data Explorer Native Files`_.

.. _`B.2 Data Explorer Native Files`:
   https://web.archive.org/web/20080808140524/http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm
.. http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF

Classes and functions
---------------------

"""
from __future__ import with_statement

import numpy
import re
from six import next
from six.moves import range


class DXclass(object):
    """'class' object as defined by OpenDX"""
    def __init__(self,classid):
        """id is the object number"""
        self.id = classid  # serial number of the object
        self.name = None   # name of the DXclass
        self.component = None   # component type
        self.D = None      # dimensions
    def write(self,file,optstring="",quote=False):
        """write the 'object' line; additional args are packed in string"""
        classid = str(self.id)
        if quote: classid = '"'+classid+'"'
        # Only use a *single* space between tokens; both chimera's and pymol's DX parser
        # does not properly implement the OpenDX specs and produces garbage with multiple
        # spaces. (Chimera 1.4.1, PyMOL 1.3)
        file.write('object '+classid+' class '+str(self.name)+' '+\
                   optstring+'\n')

    def read(self,file):
        raise NotImplementedError('Reading is currently not supported.')

    def ndformat(self,s):
        """Returns a string with as many repetitions of s as self
        has dimensions (derived from shape)"""
        return s * len(self.shape)

    def __repr__(self):
        return '<OpenDX.'+str(self.name)+' object, id='+str(self.id)+'>'


class gridpositions(DXclass):
    """OpenDX gridpositions class.

    shape     D-tuplet describing size in each dimension
    origin    coordinates of the centre of the grid cell with index 0,0,...,0
    delta     DxD array describing the deltas
    """
    def __init__(self,classid,shape=None,origin=None,delta=None,**kwargs):
        if shape is None or origin is None or delta is None:
            raise ValueError('all keyword arguments are required')
        self.id = classid
        self.name = 'gridpositions'
        self.component = 'positions'
        self.shape = numpy.asarray(shape)      # D dimensional shape
        self.origin = numpy.asarray(origin)    # D vector
        self.rank = len(self.shape)            # D === rank

        self.delta = numpy.asarray(delta)      # DxD array of grid spacings
        # gridDataFormats  actually provides a simple 1D array with the deltas because only
        # regular grids are used but the following is a reminder that OpenDX should be able
        # to handle more complicated volume elements
        if len(self.delta.shape) == 1:
            self.delta = numpy.diag(delta)
        if self.delta.shape != (self.rank, self.rank):
            # check OpenDX specs for irreg spacing if we want to implement
            # anything more complicated
            raise NotImplementedError('Only regularly spaced grids allowed, '
                                      'not delta={}'.format(self.delta))
    def write(self,file):
        DXclass.write(self,file,
                      ('counts '+self.ndformat(' %d')) % tuple(self.shape))
        file.write('origin %f %f %f\n' % tuple(self.origin))
        for delta in self.delta:
            file.write(('delta '+self.ndformat(' %f')+'\n') % tuple(delta))
    def edges(self):
        """Edges of the grid cells, origin at centre of 0,0,..,0 grid cell.

        Only works for regular, orthonormal grids.
        """
        return [self.delta[d,d] * numpy.arange(self.shape[d]+1) + self.origin[d]\
                - 0.5*self.delta[d,d]     for d in range(self.rank)]


class gridconnections(DXclass):
    """OpenDX gridconnections class"""
    def __init__(self,classid,shape=None,**kwargs):
        if shape is None:
            raise ValueError('all keyword arguments are required')
        self.id = classid
        self.name = 'gridconnections'
        self.component = 'connections'
        self.shape = numpy.asarray(shape)      # D dimensional shape
    def write(self,file):
        DXclass.write(self,file,
                      ('counts '+self.ndformat(' %d')) % tuple(self.shape))

class array(DXclass):
    """OpenDX array class"""
    def __init__(self,classid,array=None,**kwargs):
        if array is None:
            raise ValueError('array keyword argument is required')
        self.id = classid
        self.name = 'array'
        self.component = 'data'
        self.array = numpy.asarray(array)
    def write(self,file):
        DXclass.write(self,file,
                      'type float rank 0 items %d data follows' % \
                      self.array.size)
        # grid data, serialized as a C array (z fastest varying)
        # (flat iterator is equivalent to: for x: for y: for z: grid[x,y,z])
        # VMD's DX reader requires exactly 3 values per line
        values_per_line = 3
        values = self.array.flat
        while 1:
            try:
                for i in range(values_per_line):
                    file.write(str(next(values)) + "\t")
                file.write('\n')
            except StopIteration:
                file.write('\n')
                break
        file.write('attribute "dep" string "positions"\n')

class field(DXclass):
    """OpenDX container class

    The *field* is the top-level object and represents the whole
    OpenDX file. It contains a number of other objects.

    Methods overview:

    :meth:`add`
        add a component to the field
    :meth:`add_comments`
        add comments
    :meth:`write`
       write OpenDX file to file descriptor; only simple regular
       arrays are supported. File should be readable by VMD.
    :meth:`read`
        construct the field from a dx file

    Instantiated a DX object from this class and add subclasses with
    :meth:`append`.
    """
    # perhaps this should not derive from DXclass as those are
    # objects in field but a field cannot contain itself
    def __init__(self,classid='0',components=None,comments=None):
        """OpenDX object, which is build from a list of components.::

          dx = OpenDX.field('density',[gridpoints,gridconnections,array])

        :Arguments:

           *id*
               arbitrary string
           *components*
               dictionary of DXclass instances (no sanity check on the
               individual ids!) which correspond to
                  * positions
                  * connections
                  * data
           *comments*
               list of strings; each string becomes a comment line
               prefixed with '#'. Avoid newlines.

        A field must have at least the components 'positions',
        'connections', and 'data'. Those components are associated
        with objects belonging to the field. When writing a dx file
        from the field, only the required objects are dumped to the file.

        (For a more general class that can use field:
        Because there could be more objects than components, we keep a
        separate object list. When dumping the dx file, first all
        objects are written and then the field object describes its
        components. Objects are referenced by their unique id.)

        .. Warning:: uniqueness of the *id* is not checked.

        """
        if components is None:
            components = dict(positions=None,connections=None,data=None)
        if comments is None:
            comments = ['OpenDX written by gridData.OpenDX',
                        'from http://github.com/orbeckst/GridDataFormats']
        elif type(comments) is not list:
            comments = [str(comments)]
        self.id = classid       # can be an arbitrary string
        self.name = 'field'
        self.component = None   # cannot be a component of a field
        self.components = components
        self.comments= comments

    def write(self,filename):
        """Write the complete dx object to the file.

        write(filename)

        This is the simple OpenDX format which includes the data into
        the header via the 'object array ... data follows' statement.

        The format should be compatible with VMD's dx reader plugin.
        """
        # comments (VMD chokes on lines of len > 80, so truncate)
        maxcol = 80
        with open(filename,'w') as outfile:
            for line in self.comments:
                comment = '# '+str(line)
                outfile.write(comment[:maxcol]+'\n')
            # each individual object
            for component,object in self.sorted_components():
                object.write(outfile)
            # the field object itself
            DXclass.write(self,outfile,quote=True)
            for component,object in self.sorted_components():
                outfile.write('component "%s" value %s\n' % (component,str(object.id)))

    def read(self,file):
        """Read DX field from file.

        dx = OpenDX.field.read(dxfile)

        The classid is discarded and replaced with the one from the file.
        """
        DXfield = self
        p = DXParser(file)
        p.parse(DXfield)

    def add(self,component,DXobj):
        self[component] = DXobj

    def add_comment(self,comment):
        self.comments.append(comment)

    def sorted_components(self):
        """iterator that returns (component,object) in id order"""
        for component, object in \
                sorted(self.components.items(),
                       key=lambda comp_obj: comp_obj[1].id):
            yield component, object

    def histogramdd(self):
        """Return array data as (edges,grid), i.e. a numpy nD histogram."""
        shape = self.components['positions'].shape
        edges = self.components['positions'].edges()
        hist = self.components['data'].array.reshape(shape)
        return (hist,edges)

    def __getitem__(self,key):
        return self.components[key]

    def __setitem__(self,key,value):
        self.components[key] = value

    def __repr__(self):
        return '<OpenDX.field object, id='+str(self.id)+', with '+\
               str(len(self.components))+' components and '+\
               str(len(self.components))+' objects>'


#------------------------------------------------------------
# DX file parsing
#------------------------------------------------------------

class DXParseError(Exception):
    """general exception for parsing errors in DX files"""
    pass
class DXParserNoTokens(DXParseError):
    """raised when the token buffer is exhausted"""
    pass

class Token:
    # token categories (values of dx_regex must match up with these categories)
    category = {'COMMENT': ['COMMENT'],
                'WORD': ['WORD'],
                'STRING': ['QUOTEDSTRING','BARESTRING','STRING'],
                'WHITESPACE': ['WHITESPACE'],
                'INTEGER': ['INTEGER'],
                'REAL': ['REAL'],
                'NUMBER': ['INTEGER','REAL']}
    # cast functions
    cast = {'COMMENT': lambda s:re.sub(r'#\s*','',s),
            'WORD': str,
            'STRING': str, 'QUOTEDSTRING': str, 'BARESTRING': str,
            'WHITESPACE': None,
            'NUMBER': float, 'INTEGER': int, 'REAL': float}

    def __init__(self,code,text):
        self.code = code    # store raw code
        self.text = text
    def equals(self,v):
        return self.text == v
    def iscode(self,code):
        return self.code in self.category[code]  # use many -> 1 mappings
    def value(self,ascode=None):
        """Return text cast to the correct type or the selected type"""
        if ascode is None:
            ascode = self.code
        return self.cast[ascode](self.text)
    def __repr__(self):
        return '<token '+str(self.code)+','+str(self.value())+'>'

class DXInitObject(object):
    """Storage class that holds data to initialize one of the 'real'
    classes such as OpenDX.array, OpenDX.gridconnections, ...

    All variables are stored in args which will be turned into the
    arguments for the DX class.
    """
    DXclasses = {'gridpositions':gridpositions,
                 'gridconnections':gridconnections,
                 'array':array, 'field':field,
                 }

    def __init__(self,classtype,classid):
        self.type = classtype
        self.id = classid
        self.args = dict()
    def initialize(self):
        """Initialize the corresponding DXclass from the data.

        class = DXInitObject.initialize()
        """
        return self.DXclasses[self.type](self.id,**self.args)
    def __getitem__(self,k):
        return self.args[k]
    def __setitem__(self,k,v):
        self.args[k] = v
    def __repr__(self):
        return '<DXInitObject instance type='+str(self.type)+', id='+str(self.id)+'>'

class DXParser(object):
    """Brain-dead baroque implementation to read a simple (VMD) dx file.

    Requires a OpenDX.field instance.

    1) scan for 'object' lines:
       'object' id 'class' class  [data]
       [data ...]
    2) parse data according to class
    3) construct dx field from classes
    """

    # the regexes must match with the categories defined in the Token class
    dx_regex = re.compile(r"""
    (?P<COMMENT>\#.*$)            # comment (until end of line)
    |(?P<WORD>(object|class|counts|origin|delta|type|counts|rank|items|data))
    |"(?P<QUOTEDSTRING>[^\"]*)"   # string in double quotes  (quotes removed)
    |(?P<WHITESPACE>\s+)          # white space
    |(?P<REAL>[-+]?               # true real number (decimal point or
          (\d+\.\d*([eE][-+]\d+)?)  # scientific notation)
          |(\d*\.\d+([eE][-+]\d+)?)
          |(\d[eE][-+]\d+))
    |(?P<INTEGER>[-+]?\d+)       # integer
    |(?P<BARESTRING>[a-zA-Z_][^\s\#\"]+) # unquoted strings, starting with non-numeric
    """, re.VERBOSE)


    def __init__(self,filename):
        """Setup a parser for a simple DX file (from VMD)

        >>> DXfield_object = OpenDX.field(id)
        >>> p = DXparser('bulk.dx')
        >>> p.parse(DXfield_object)

        The field object will be completely rewritten (including the
        id if one is found in the input file. The input files
        component layout is currently ignored.

        Note that quotes are removed from quoted strings.
        """
        self.filename = filename
        self.field = field('grid data',comments=['filename: '+self.filename])
        # other variables are initialised every time parse() is called

        self.parsers = {'general':self.__general,
                        'comment':self.__comment, 'object':self.__object,
                        'gridpositions':self.__gridpositions,
                        'gridconnections':self.__gridconnections,
                        'array':self.__array, 'field':self.__field,
                        }


    def parse(self,DXfield):
        """Parse the dx file and construct a DX field object with component classes.

        A :class:`field` instance *DXfield* must be provided to be
        filled by the parser::

           DXfield_object = OpenDX.field(*args)
           parse(DXfield_object)

        A tokenizer turns the dx file into a stream of tokens. A
        hierarchy of parsers examines the stream. The level-0 parser
        ('general') distinguishes comments and objects (level-1). The
        object parser calls level-3 parsers depending on the object
        found. The basic idea is that of a 'state machine'. There is
        one parser active at any time. The main loop is the general
        parser.

        * Constructing the dx objects with classtype and classid is
          not implemented yet.
        * Unknown tokens raise an exception.
        """

        self.DXfield = DXfield              # OpenDX.field (used by comment parser)
        self.currentobject = None           # containers for data
        self.objects = []                   # |
        self.tokens = []                    # token buffer
        with open(self.filename,'r') as self.dxfile:
            self.use_parser('general')      # parse the whole file and populate self.objects

        # assemble field from objects
        for o in self.objects:
            if o.type == 'field':
                # Almost ignore the field object; VMD, for instance,
                # does not write components. To make this work
                # seamlessly I have to think harder how to organize
                # and use the data, eg preping the field object
                # properly and the initializing. Probably should also
                # check uniqueness of ids etc.
                DXfield.id = o.id
                continue
            c = o.initialize()
            self.DXfield.add(c.component,c)

        # free space
        del self.currentobject, self.objects



    def __general(self):
        """Level-0 parser and main loop.

        Look for a token that matches a level-1 parser and hand over control."""
        while 1:                            # main loop
            try:
                tok = self.__peek()         # only peek, apply_parser() will consume
            except DXParserNoTokens:
                # save previous DXInitObject
                # (kludge in here as the last level-2 parser usually does not return
                # via the object parser)
                if self.currentobject and self.currentobject not in self.objects:
                    self.objects.append(self.currentobject)
                return                      # stop parsing and finish
            # decision branches for all level-1 parsers:
            # (the only way to get out of the lower level parsers!)
            if tok.iscode('COMMENT'):
                self.set_parser('comment')  # switch the state
            elif tok.iscode('WORD') and tok.equals('object'):
                self.set_parser('object')   # switch the state
            elif self.__parser is self.__general:
                # Either a level-2 parser screwed up or some level-1
                # construct is not implemented.  (Note: this elif can
                # be only reached at the beginning or after comments;
                # later we never formally switch back to __general
                # (would create inifinite loop)
                raise DXParseError('Unknown level-1 construct at '+str(tok))

            self.apply_parser()     # hand over to new parser
                                    # (possibly been set further down the hierarchy!)

    # Level-1 parser
    def __comment(self):
        """Level-1 parser for comments.

        pattern: #.*
        Append comment (with initial '# ' stripped) to all comments.
        """
        tok = self.__consume()
        self.DXfield.add_comment(tok.value())
        self.set_parser('general')   # switch back to general parser

    def __object(self):
        """Level-1 parser for objects.

        pattern: 'object' id 'class' type ...

        id ::=   integer|string|'"'white space string'"'
        type ::= string
        """
        self.__consume()                    # 'object'
        classid = self.__consume().text
        word = self.__consume().text
        if word != "class":
            raise DXParseError("reserved word %s should have been 'class'." % word)
        # save previous DXInitObject
        if self.currentobject:
            self.objects.append(self.currentobject)
        # setup new DXInitObject
        classtype = self.__consume().text
        self.currentobject = DXInitObject(classtype=classtype,classid=classid)

        self.use_parser(classtype)

    # Level-2 parser (object parsers)
    def __gridpositions(self):
        """Level-2 parser for gridpositions.

        pattern:
        object 1 class gridpositions counts 97 93 99
        origin -46.5 -45.5 -48.5
        delta 1 0 0
        delta 0 1 0
        delta 0 0 1
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('counts'):
            shape = []
            try:
                while self.__peek().iscode('INTEGER'):
                    tok = self.__consume()
                    shape.append(tok.value())
            except DXParserNoTokens:
                pass
            if len(shape) == 0:
                raise DXParseError('gridpositions: no shape parameters')
            self.currentobject['shape'] = shape
        elif tok.equals('origin'):
            origin = []
            try:
                while (self.__peek().iscode('INTEGER') or
                       self.__peek().iscode('REAL')):
                    tok = self.__consume()
                    origin.append(tok.value())
            except DXParserNoTokens:
                pass
            if len(origin) == 0:
                raise DXParseError('gridpositions: no origin parameters')
            self.currentobject['origin'] = origin
        elif tok.equals('delta'):
            d = []
            try:
                while (self.__peek().iscode('INTEGER') or
                       self.__peek().iscode('REAL')):
                    tok = self.__consume()
                    d.append(tok.value())
            except DXParserNoTokens:
                pass
            if len(d) == 0:
                raise DXParseError('gridpositions: missing delta parameters')
            try:
                self.currentobject['delta'].append(d)
            except KeyError:
                self.currentobject['delta'] = [d]
        else:
            raise DXParseError('gridpositions: '+str(tok)+' not recognized.')


    def __gridconnections(self):
        """Level-2 parser for gridconnections.

        pattern:
        object 2 class gridconnections counts 97 93 99
        """

        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('counts'):
            shape = []
            try:
                while self.__peek().iscode('INTEGER'):
                    tok = self.__consume()
                    shape.append(tok.value())
            except DXParserNoTokens:
                pass
            if len(shape) == 0:
                raise DXParseError('gridconnections: no shape parameters')
            self.currentobject['shape'] = shape
        else:
            raise DXParseError('gridconnections: '+str(tok)+' not recognized.')


    def __array(self):
        """Level-2 parser for arrays.

        pattern:
        object 3 class array type double rank 0 items 12 data follows
        0 2 0
        0 0 3.6
        0 -2.0 1e-12
        +4.534e+01 .34534 0.43654
        attribute "dep" string "positions"
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('type'):
            tok = self.__consume()
            if not tok.iscode('STRING'):
                raise DXParseError('array: type was "%s", not a string.'%\
                                   tok.text)
            self.currentobject['type'] = tok.value()
        elif tok.equals('rank'):
            tok = self.__consume()
            if not tok.iscode('INTEGER'):
                raise DXParseError('array: rank was "%s", not an integer.'%\
                                   tok.text)
            self.currentobject['rank'] = tok.value()
        elif tok.equals('items'):
            tok = self.__consume()
            if not tok.iscode('INTEGER'):
                raise DXParseError('array: items was "%s", not an integer.'%\
                                   tok.text)
            self.currentobject['size'] = tok.value()
        elif tok.equals('data'):
            tok = self.__consume()
            if not tok.iscode('STRING'):
                raise DXParseError('array: data was "%s", not a string.'%\
                                   tok.text)
            if tok.text != 'follows':
                raise NotImplementedError(\
                            'array: Only the "data follows header" format is supported.')
            if not self.currentobject['size']:
                raise DXParseError("array: missing number of items")
            self.currentobject['array'] = [self.__consume().value('REAL') \
                                           for i in range(self.currentobject['size'])]
        elif tok.equals('attribute'):
            # not used at the moment
            attribute = self.__consume().value()
            if not self.__consume().equals('string'):
                raise DXParseError('array: "string" expected.')
            value = self.__consume().value()
        else:
            raise DXParseError('array: '+str(tok)+' not recognized.')

    def __field(self):
        """Level-2 parser for a DX field object.

        pattern:
        object "site map 1" class field
        component "positions" value 1
        component "connections" value 2
        component "data" value 3
        """
        try:
            tok = self.__consume()
        except DXParserNoTokens:
            return

        if tok.equals('component'):
            component = self.__consume().value()
            if not self.__consume().equals('value'):
                raise DXParseError('field: "value" expected')
            classid = self.__consume().value()
            try:
                self.currentobject['components'][component] = classid
            except KeyError:
                self.currentobject['components'] = {component:classid}
        else:
            raise DXParseError('field: '+str(tok)+' not recognized.')

    # parser routines independent of the dx classes
    # (with ideas from MDAnalysis.Selection and
    # http://effbot.org/zone/xml-scanner.htm)

    def use_parser(self,parsername):
        """Set parsername as the current parser and apply it."""
        self.__parser = self.parsers[parsername]
        self.__parser()
    def set_parser(self,parsername):
        """Set parsername as the current parser."""
        self.__parser = self.parsers[parsername]
    def apply_parser(self):
        """Apply the current parser to the token stream."""
        self.__parser()

    def __tokenize(self,string):
        """Split s into tokens and update the token buffer.

        __tokenize(string)

        New tokens are appended to the token buffer, discarding white
        space.  Based on http://effbot.org/zone/xml-scanner.htm
        """
        for m in self.dx_regex.finditer(string.strip()):
            code = m.lastgroup
            text = m.group(m.lastgroup)
            tok = Token(code,text)
            if not tok.iscode('WHITESPACE'):
                 self.tokens.append(tok)
                 # print "DEBUG tokenize: "+str(tok)

    def __refill_tokenbuffer(self):
        """Add a new tokenized line from the file to the token buffer.

        __refill_tokenbuffer()

        Only reads a new line if the buffer is empty. It is safe to
        call it repeatedly.

        At end of file, method returns empty strings and it is up to
        __peek and __consume to flag the end of the stream.
        """
        if len(self.tokens) == 0:
            self.__tokenize(self.dxfile.readline())

    def __peek(self):
        self.__refill_tokenbuffer()
        try:
            return self.tokens[0]
        except IndexError:
            raise DXParserNoTokens

    def __consume(self,):
        """Get the next token from the buffer and remove it/them.

        try:
          while 1:
             token = __consume()
        except DXParserNoTokens:
          pass
        """
        self.__refill_tokenbuffer()
        #print "DEBUG consume: "+str(self.__parser)+' '+str(self.__peek())
        try:
            return self.tokens.pop(0)  # singlet
        except IndexError:
            raise DXParserNoTokens
