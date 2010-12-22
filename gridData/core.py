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
import OpenDX, gOpenMol

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

    def __init__(self,grid=None,edges=None,origin=None,delta=None, metadata=None, **kwargs):
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
          interpolation_spline_order
            order of interpolation function for resampling; cubic splines = 3 [3]
        """
        # file formats are guess from extension == lower case key
        self._exporters = {'DX': self._export_dx,
                           'PICKLE': self._export_python,
                           'PYTHON': self._export_python,  # compatibility
                           }
        self._loaders = {'DX': self._load_dx,
                         'PLT': self._load_plt,
                         'PICKLE': self._load_python,
                         'PYTHON': self._load_python,      # compatibility
                         }

        if metadata is None: 
            metadata = {}
        self.metadata = metadata     # use this to record arbitrary data
        self.__interpolated = None   # cache for interpolated grid
        self.__interpolation_spline_order = kwargs.pop('interpolation_spline_order', 3)
        self.interpolation_cval = None       # default to using min(grid)

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
            self.grid = numpy.asarray(grid)
            self._update()
        else:
            # empty, must manually populate with load()
            #print "Setting up empty grid object. Use Grid.load(filename)."
            pass

    def interpolation_spline_order():
        """Order of the B-spline interpolation of the data.

        3 = cubic; 4 & 5 are also supported

        Only choose values that are acceptable to :func:`scipy.ndimage.spline_filter`!
        """
        def fget(self):
            return self.__interpolation_spline_order
        def fset(self, x):
            # As we cache the interpolation function, we need to rebuild the cache
            # whenever the interpolation order changes: this is handled by _update()
            self.__interpolation_spline_order = x
            self._update()
        return locals()
    interpolation_spline_order = property(**interpolation_spline_order())

    def resample(self, edges):
        """Resample data to a new grid with edges *edges*.

          resample(edges) --> Grid

        or 

          resample(otherGrid) --> Grid

        The order of the interpolation is set by
        :attr:`Grid.interpolation_spline_order`.
        """
        try:
            edges = edges.edges  # can also supply another Grid
        except AttributeError:
            pass
        midpoints = self._midpoints(edges)
        coordinates = ndmeshgrid(*midpoints)
        newgrid = self.interpolated(*coordinates)  # feed a meshgrid to generate all points
        return Grid(newgrid, edges)

    def resample_factor(self, factor):
        """Resample to a new regular grid with factor*oldN cells along each dimension."""
        from itertools import izip
        # new number of edges N' = (N-1)*f + 1
        newlengths = [(N-1)*float(factor) + 1 for N in self._len_edges()]
        edges = [numpy.linspace(start,stop,num=N,endpoint=True) for (start,stop,N) in 
                 izip(self._min_edges(), self._max_edges(), newlengths)]
        return self.resample(edges)

    def _edgify(self, midpoints):
        """Return edges, given midpoints."""
        m = numpy.asarray(midpoints)
        return numpy.concatenate([[m[0] - 0.5*(m[1]-m[0])], m, [m[-1] + 0.5*(m[-1]-m[-2])]])

    def _update(self):
        """compute/update all derived data

        Grid._update()

        Can be called without harm and is idem-potent.

        Updates these attributes and methods:
           :attr:`origin`
              the center of the cell with index 0,0,0
           :attr:`midpoints`
              centre coordinate of each grid cell
           :meth:`interpolated`
              spline interpolation function that can generated a value for
              coordinate
        """
        self.delta = numpy.diag(
            map(lambda e: (e[-1] - e[0])/(len(e)-1), self.edges) )
        self.midpoints = self._midpoints(self.edges)
        self.origin = numpy.array(map(lambda m: m[0], self.midpoints))
        if not self.__interpolated is None:
            # only update if we are using it
            self.__interpolated = self._interpolationFunctionFactory()

    @property
    def interpolated(self):
        """B-spline function over the data grid(x,y,z).

           interpolated([x1,x2,...],[y1,y2,...],[z1,z2,...]) -> F[x1,y1,z1],F[x2,y2,z2],...

        The interpolation order is set in :attr:`Grid.interpolation_spline_order`.

        The interpolated function is computed once and is cached for better
        performance. Whenever :attr:`~Grid.interpolation_spline_order` is
        modified, :meth:`Grid.interpolated` is recomputed.

        The value for unknown data is set in :attr:`Grid.interpolation_cval`
        (TODO: also recompute when interpolation_cval value is changed.)

        Example usage for resampling::
           >>> XX,YY,ZZ = numpy.mgrid[40:75:0.5, 96:150:0.5, 20:50:0.5]
           >>> FF = interpolated(XX,YY,ZZ)            
        """
        if self.__interpolated is None:
            self.__interpolated = self._interpolationFunctionFactory()
        return self.__interpolated        

    def _map_edges(self, func, edges=None):
        if edges is None:
            edges = self.edges
        return [func(e) for e in edges]

    def _midpoints(self, edges=None):
        return self._map_edges(lambda e: 0.5*(e[:-1] + e[1:]), edges=edges)

    def _len_edges(self, edges=None):
        return self._map_edges(len, edges=edges)    

    def _min_edges(self, edges=None):
        return self._map_edges(numpy.min, edges=edges)

    def _max_edges(self, edges=None):
        return self._map_edges(numpy.max, edges=edges)


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

    def _load_dx(self, filename):
        """Initializes Grid from a OpenDX file."""
        
        dx = OpenDX.field(0)
        dx.read(filename)
        grid,edges = dx.histogramdd()
        self.__init__(grid=grid,edges=edges,metadata=self.metadata)
    
    def _load_plt(self, filename):
        """Initialize Grid from gOpenMol plt file."""
        g = gOpenMol.Plt()
        g.read(filename)
        grid,edges = g.histogramdd()
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
                            "It must be a scalar or a grid with identical edges. "
                            "Use Grid.resample(other.edges) to make a new grid that is "
                            "compatible with other.")
        return True

    def _interpolationFunctionFactory(self,spline_order=None,cval=None):
        """Returns a function F(x,y,z) that interpolates any values on the grid.

        _interpolationFunctionFactory(self,spline_order=3,cval=None) --> F

        *cval* is set to :meth:`Grid.grid.min`. *cval* cannot be chosen too
        large or too small or NaN because otherwise the spline interpolation
        breaks down near that region and produces wild oscillations.

        .. Note:: Only correct for equally spaced values (i.e. regular edges with
                  constant delta).
        .. SeeAlso:: http://www.scipy.org/Cookbook/Interpolation
        """
        # for scipy >=0.9: should use scipy.interpolate.griddata
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html#scipy.interpolate.griddata
        # (does it work for nD?)

        from scipy import ndimage

        if spline_order is None:
            # must be compatible with whatever :func:`scipy.ndimage.spline_filter` takes.
            spline_order = self.interpolation_spline_order
        if cval is None:
            cval = self.interpolation_cval

        data = self.grid
        if cval is None:
            cval = data.min()
        try:
            # masked arrays
            _data = data.filled(cval)   # fill with min; hopefully keeps spline happy
        except AttributeError:
            _data = data

        coeffs = ndimage.spline_filter(_data,order=spline_order)
        x0 = self.origin
        dx = self.delta.diagonal()    # fixed dx required!!
        def _transform(cnew, c0, dc):
            return (numpy.atleast_1d(cnew) - c0)/dc
        def interpolatedF(*coordinates):
            """B-spline function over the data grid(x,y,z).

            interpolatedF([x1,x2,...],[y1,y2,...],[z1,z2,...]) -> F[x1,y1,z1],F[x2,y2,z2],...

            Example usage for resampling::
              >>> XX,YY,ZZ = numpy.mgrid[40:75:0.5, 96:150:0.5, 20:50:0.5]
              >>> FF = _interpolationFunction(XX,YY,ZZ)            
            """
            _coordinates = numpy.array(
                [_transform(coordinates[i], x0[i], dx[i]) for i in xrange(len(coordinates))])
            return ndimage.map_coordinates(coeffs, _coordinates, prefilter=False, 
                                           mode='nearest',cval=cval)
        # mode='wrap' would be ideal but is broken: http://projects.scipy.org/scipy/ticket/796
        return interpolatedF            
                

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
        try:
            bins = self.grid.shape
        except AttributeError:
            bins = "no"
        return '<Grid with '+str(bins)+' bins>'

def ndmeshgrid(*arrs):
    """Return a mesh grid for N dimensions.

    The input are N arrays, each of which contains the values along one axis of
    the coordinate system. The arrays do not have to have the same number of
    entries. The function returns arrays that can be fed into numpy functions
    so that they produce values for *all* points spanned by the axes *arrs*.

    Original from 
    http://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d and fixed.

    .. SeeAlso: :func:`numpy.meshgrid` for the 2D case.
    """
    #arrs = tuple(reversed(arrs)) <-- wrong on stackoverflow.com
    arrs = tuple(arrs)
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []    
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = numpy.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j) 
        ans.append(arr2)

    return tuple(ans)
