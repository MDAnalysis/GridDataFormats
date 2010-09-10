# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser Public License, version 3 or later.

"""
:mod:`gridData.core` --- Core functionality for storing n-D grids
=================================================================
 
Classes and functions that are independent of the grid data
format. In particular this module contains the :class:`Grid` class that acts as
a universal constructor for specific formats::

 g = Grid(**kwargs)           # construct
 g.export(filename, format)   # export to the desired format

Some formats can also be read::

 g = Grid()                   # make an empty Grid
 g.load(filename)             # populate with data from filename

"""

import os
import warnings
import cPickle
import numpy
import OpenDX

from gridData import gridDataWarning

def _grid(x):
    """Access the underlying ndarray of a Grid object or return the object itself"""
    try:
        return x.grid
    except AttributeError:
        return x

class Grid(object):
    """Class to manage a multidimensional grid object.

    The export(format='dx') method always exports a 3D object, the
    rest should work for an array of any dimension.

    The grid (Grid.grid) can be manipulated as a standard numpy
    array. 

    The attribute Grid.metadata holds a user-defined dictionary that
    can be used to annotate the data. It is saved with save().
    """
    default_format = 'DX'

    def __init__(self,grid=None,edges=None,origin=None,delta=None, metadata=None):
        """
        Create a Grid object from data.

        From a numpy.histogramdd()::
          grid,edges = numpy.histogramdd(...)
          g = Grid(grid,edges=edges)

        From an arbitrary grid::
          g = Grid(grid,origin=origin,delta=delta)

        From a saved file::
          g = Grid(filename)
        or
          g = Grid()
          g.load(filename)
 
        :Arguments:
          grid       
            histogram or density, defined on numpy nD array
          edges
            list of arrays, the lower and upper bin edges along the axes
            (both are output by numpy.histogramdd())
          origin
            cartesian coordinates of the center of grid[0,0,...,0]
          delta
            Either n x n array containing the cell lengths in each dimension, 
            or n x 1 array for rectangular arrays.
          metadata
            a user defined dictionary of arbitrary values
            associated with the density; the class does not touch
            metadata[] but stores it with save()
        """
        self._exporters = {'DX': self._export_dx,
                           'PICKLE': self._export_python,
                           'PYTHON': self._export_python,  # compatibility
                           }
        self._loaders = {'DX': self._load_dx,
                         'PICKLE': self._load_python,
                         'PYTHON': self._load_python,      # compatibility
                         }


        if metadata is None: 
            metadata = {}
        self.metadata = metadata     # use this to record arbitrary data

        if type(grid) is str:
            # read from a file
            self.load(grid)
        elif not (grid is None or edges is None):
            # set up from histogramdd-type data
            self.grid = numpy.asarray(grid)
            self.edges = edges
            self._update()
        elif not (grid is None or origin is None or delta is None):
            # setup from generic data
            origin = numpy.squeeze(origin)
            delta = numpy.squeeze(delta)
            N = grid.ndim
            assert(N == len(origin))
            if delta.shape == (N,N):
                if numpy.any(delta - numpy.diag(delta)):
                    raise NotImplementedError("Non-rectangular grids are not supported.")
            elif delta.shape == (N,):
                delta = numpy.diagflat(delta)
            elif delta.shape == ():
                delta = numpy.diagflat(N*[delta])
            else:
                raise ValueError('delta = %r has the wrong shape' % delta)
            # note that origin is CENTER so edges must be shifted by -0.5*delta
            self.edges = [origin[dim] + (numpy.arange(m+1) - 0.5) * delta[dim,dim] 
                          for dim,m in enumerate(grid.shape)]
            self.grid = grid
            self._update()
        else:
            # empty, must manually populate with load()
            #print "Setting up empty grid object. Use Grid.load(filename)."
            pass

    def _update(self):
        """compute/update all derived data

        Grid._update()

        Can be called without harm and is idem-potent.

        origin  is the center of the cell with index 0,0,0
        """
        self.delta = numpy.diag(
            map(lambda e: (e[-1] - e[0])/(len(e)-1), self.edges) )
        self.midpoints = map(lambda e: 0.5 * (e[:-1] + e[1:]), self.edges)
        self.origin = map(lambda m: m[0], self.midpoints)

    def _guess_format(self, filename, format=None, export=True):
        if export:
            available = self._exporters
        else:
            available = self._loaders
        if format is None:
            format = os.path.splitext(filename)[1][1:]
        format = format.upper()
        if not format:
            format = self.default_format
        if not format in available:
            raise ValueError("File format %r not available, choose one of %r"\
                                 % (format, available.keys()))
        return format

    def _get_exporter(self, filename, format=None):
        return self._exporters[self._guess_format(filename, format=format, export=True)]

    def _get_loader(self, filename, format=None):
        return self._loaders[self._guess_format(filename, format=format, export=False)]

    def load(self,filename, format=None):
        """Load saved (pickled or dx) grid and edges from <filename>.pickle

           Grid.load(<filename>.pickle)
           Grid.load(<filename>.dx)

        The load() method calls the class's constructor method and
        completely resets all values, based on the loaded data.
        """
        loader = self._get_loader(filename, format=format)
        loader(filename)

    def _load_python(self,filename):
        f = open(filename,'rb')
        try:
            saved = cPickle.load(f)
        finally:
            f.close()
        self.__init__(grid=saved['grid'],edges=saved['edges'],metadata=saved['metadata'])
        del saved

    def _load_dx(self,dxfile):
        """Initializes Grid from a OpenDX file."""
        
        dx = OpenDX.field(0)
        dx.read(dxfile)
        grid,edges = dx.histogramdd()
        self.__init__(grid=grid,edges=edges,metadata=self.metadata)
    
    def export(self,filename,format=None):
        """export density to file using the given format; use 'dx' for visualization.

        export(filename=<filename>,format=<format>)

        The format can also be deduced from the suffix of the filename
        though the *format* keyword takes precedence.

        The default format for export() is 'dx'.
        
        Only implemented formats:

        dx        OpenDX
        pickle    pickle (use Grid.load(filename) to restore); Grid.save()
                  is simpler than export(format='python').
        """
        exporter = self._get_exporter(filename, format=format)
        exporter(filename)

    def _export_python(self,filename):
        """Pickle the Grid object

        The object is dumped as a dictionary with grid and edges: This
        is sufficient to recreate the grid object with __init__().
        """
        root, ext = os.path.splitext(filename)
        filename = root + ".pickle"
        
        data = dict(grid=self.grid,edges=self.edges,metadata=self.metadata)
        f = open(filename,'wb')
        try:
            cPickle.dump(data,f,cPickle.HIGHEST_PROTOCOL)
        finally:
            f.close()
        del data

    def _export_dx(self,filename):
        """Export the density grid to an OpenDX file. The file format
        is the simplest regular grid array and it is also understood
        by VMD's and PyMOL's DX reader.

        For the file format see
        http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF
        """
        root, ext = os.path.splitext(filename)
        filename = root + '.dx'

        comments = [
            'OpenDX density file written by gridDataFormats.Grid.export()',
            'File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF',
            'Data are embedded in the header and tied to the grid positions.',
            'Data is written in C array order: In grid[x,y,z] the axis z is fastest',
            'varying, then y, then finally x, i.e. z is the innermost loop.']

        # write metadata in comments section
        if self.metadata:
            comments.append('Meta data stored with the python Grid object:')
        for k in self.metadata:
            comments.append('   '+str(k)+' = '+str(self.metadata[k]))
        comments.append('(Note: the VMD dx-reader chokes on comments below this line)')

        components = dict(
            positions = OpenDX.gridpositions(1,self.grid.shape,self.origin,self.delta),
            connections = OpenDX.gridconnections(2,self.grid.shape),
            data = OpenDX.array(3,self.grid),
            )
        dx = OpenDX.field('density',components=components,comments=comments)
        dx.write(filename)

    def save(self,filename):
        """Save a grid object to <filename>.pickle

           Grid.save(filename)

        Internally, this calls Grid.export(filename,format="python"). A grid can be
        regenerated from the saved data with

           g = Grid(filename=<filename>)

        """
        self.export(filename,format="pickle")

    def centers(self):
        """Returns the coordinates of the centers of all grid cells as an iterator."""
        # crappy
        for idx in numpy.ndindex(self.grid.shape):
            # TODO: CHECK that this delta*(i,j,k) is really correct, even for non-diagonal delta
            # NOTE: origin is center of (0,0,0) (and already has index offset by 0.5)
            yield numpy.sum(self.delta * numpy.asarray(idx), axis=0) + self.origin

    def check_compatible(self, other):
        """Check if *other* can be used in an arithmetic operation.

        1) *other* is a scalar
        2) *other* is a grid defined on the same edges
        
        :Raises: :exc:`TypeError` if not compatible.
        """
        if not (numpy.isscalar(other) or 
                numpy.all(numpy.concatenate(self.edges) == numpy.concatenate(other.edges))):
            raise TypeError("The argument can not be arithmetically combined with the grid. "
                            "It must be a scalar or a grid with identical edges.")
        return True

    # basic arithmetic (left and right associative so that Grid1 + Grid2 but also
    # 3 * Grid and Grid/0.5 work)

    def __add__(self, other):
        """Return a new :class:`Grid` with the point-wise sum of the data.

        g.__add__(h) <==> g + h

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(self.grid + _grid(other), edges=self.edges)
    
    def __sub__(self, other):
        """Return a new :class:`Grid` with the point-wise difference of the data.

        g.__sub__(h) <==> g - h

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(self.grid - _grid(other), edges=self.edges)

    def __mul__(self, other):
        """Return a new :class:`Grid` with the point-wise product of the data.

        g.__mul__(h) <==> g * h

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(self.grid * _grid(other), edges=self.edges)

    def __div__(self, other):
        """Return a new :class:`Grid` with the point-wise quotient of the data.

        g.__div__(h) <==> g/h

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(self.grid / _grid(other), edges=self.edges)

    def __pow__(self, other):
        """Return a new :class:`Grid` with the point-wise power of the data.

        g.__pow__(h) <==> numpy.power(g, h)

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(numpy.power(self.grid, _grid(other)), edges=self.edges)

    def __radd__(self, other):
        """Return a new :class:`Grid` with the point-wise sum of the data.

        g.__add__(h) <==> h + g

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(_grid(other) + self.grid, edges=self.edges)

    def __rsub__(self, other):
        """Return a new :class:`Grid` with the point-wise difference of the data.

        g.__sub__(h) <==> h - g

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(_grid(other) - self.grid, edges=self.edges)

    def __rmul__(self, other):
        """Return a new :class:`Grid` with the point-wise product of the data.

        g.__mul__(h) <==> h * g

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(_grid(other) * self.grid, edges=self.edges)

    def __rdiv__(self, other):
        """Return a new :class:`Grid` with the point-wise quotient of the data.

        g.__div__(h) <==> h/g

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(_grid(other) / self.grid, edges=self.edges)

    def __rpow__(self, other):
        """Return a new :class:`Grid` with the point-wise power of the data.

        g.__pow__(h) <==> numpy.power(h, g)

        :Returns: :class:`Grid`
        """
        self.check_compatible(other)
        return Grid(numpy.power(_grid(other), self.grid), edges=self.edges)

    def __repr__(self):
        return '<Grid with '+str(self.grid.shape)+' bins>'
