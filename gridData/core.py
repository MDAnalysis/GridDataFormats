# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
"""
:mod:`gridData.core` --- Core functionality for storing n-D grids
=================================================================

The :mod:`core` module contains classes and functions that are
independent of the grid data format. In particular this module
contains the :class:`Grid` class that acts as a universal constructor
for specific formats::

 g = Grid(**kwargs)           # construct
 g.export(filename, format)   # export to the desired format

Some formats can also be read::

 g = Grid()                   # make an empty Grid
 g.load(filename)             # populate with data from filename


Classes and functions
---------------------

"""

import os
from six.moves import cPickle, range, zip
import numpy

# For interpolated grids: need scipy.ndimage but we import it only when needed:
# import scipy

from . import OpenDX
from . import gOpenMol
from . import CCP4


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

    def __init__(self, grid=None, edges=None, origin=None, delta=None,
                 metadata={}, interpolation_spline_order=3):
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
        self._exporters = {
            'DX': self._export_dx,
            'PICKLE': self._export_python,
            'PYTHON': self._export_python,  # compatibility
        }
        self._loaders = {
            'CCP4': self._load_cpp4,
            'DX': self._load_dx,
            'PLT': self._load_plt,
            'PICKLE': self._load_python,
            'PYTHON': self._load_python,  # compatibility
        }

        self.metadata = metadata
        self.__interpolated = None  # cache for interpolated grid
        self.__interpolation_spline_order = interpolation_spline_order
        self.interpolation_cval = None  # default to using min(grid)

        if type(grid) is str:
            self.load(grid)
        elif not (grid is None or edges is None):
            # set up from histogramdd-type data
            self.grid = numpy.asarray(grid)
            self.edges = edges
            self._update()
        elif not (grid is None or origin is None or delta is None):
            # setup from generic data
            origin = numpy.asarray(origin)
            delta = numpy.asarray(delta)
            if len(origin) != grid.ndim:
                raise TypeError(
                    "Dimension of origin is not the same as grid dimension.")
            if delta.shape == () and numpy.isreal(delta):
                delta = numpy.ones(grid.ndim) * delta
            elif delta.ndim > 1:
                raise NotImplementedError(
                    "Non-rectangular grids are not supported.")
            elif len(delta) != grid.ndim:
                raise TypeError("delta should be scalar or array-like of"
                                "len(grid.ndim)")
            # note that origin is CENTER so edges must be shifted by -0.5*delta
            self.edges = [origin[dim] +
                          (numpy.arange(m + 1) - 0.5) * delta[dim]
                          for dim, m in enumerate(grid.shape)]
            self.grid = numpy.asarray(grid)
            self._update()
        else:
            # empty, must manually populate with load()
            # print "Setting up empty grid object. Use Grid.load(filename)."
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
        # feed a meshgrid to generate all points
        newgrid = self.interpolated(*coordinates)
        return Grid(newgrid, edges)

    def resample_factor(self, factor):
        """Resample to a new regular grid with factor*oldN cells along
        each dimension."""
        # new number of edges N' = (N-1)*f + 1
        newlengths = [(N - 1) * float(factor) + 1 for N in self._len_edges()]
        edges = [numpy.linspace(start, stop, num=N, endpoint=True)
                 for (start, stop, N) in
                 zip(self._min_edges(), self._max_edges(), newlengths)]
        return self.resample(edges)

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
        self.delta = numpy.array(list(
            map(lambda e: (e[-1] - e[0]) / (len(e) - 1), self.edges)))
        self.midpoints = self._midpoints(self.edges)
        self.origin = numpy.array(list(map(lambda m: m[0], self.midpoints)))
        if self.__interpolated is not None:
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
        return self._map_edges(lambda e: 0.5 * (e[:-1] + e[1:]), edges=edges)

    def _len_edges(self, edges=None):
        return self._map_edges(len, edges=edges)

    def _min_edges(self, edges=None):
        return self._map_edges(numpy.min, edges=edges)

    def _max_edges(self, edges=None):
        return self._map_edges(numpy.max, edges=edges)

    def _guess_format(self, filename, file_format=None, export=True):
        if export:
            available = self._exporters
        else:
            available = self._loaders
        if file_format is None:
            file_format = os.path.splitext(filename)[1][1:]
        file_format = file_format.upper()
        if not file_format:
            file_format = self.default_format
        if file_format not in available:
            raise ValueError(
                "File format {} not available, choose one of {}".format(
                    file_format, available.keys()))
        return file_format

    def _get_exporter(self, filename, file_format=None):
        return self._exporters[self._guess_format(filename,
                                                  file_format=file_format,
                                                  export=True)]

    def _get_loader(self, filename, file_format=None):
        return self._loaders[self._guess_format(filename,
                                                file_format=file_format,
                                                export=False)]

    def load(self, filename, file_format=None):
        """Load saved (pickled or dx) grid and edges from <filename>.pickle

           Grid.load(<filename>.pickle)
           Grid.load(<filename>.dx)

        The load() method calls the class's constructor method and
        completely resets all values, based on the loaded data.
        """
        loader = self._get_loader(filename, file_format=file_format)
        loader(filename)

    def _load_python(self, filename):
        with open(filename, 'rb') as f:
            saved = cPickle.load(f)
        self.__init__(grid=saved['grid'],
                      edges=saved['edges'],
                      metadata=saved['metadata'])

    def _load_cpp4(self, filename):
        """Initializes Grid from a CCP4 file."""
        ccp4 = CCP4.CCP4()
        ccp4.read(filename)
        grid, edges = ccp4.histogramdd()
        self.__init__(grid=grid, edges=edges, metadata=self.metadata)

    def _load_dx(self, filename):
        """Initializes Grid from a OpenDX file."""
        dx = OpenDX.field(0)
        dx.read(filename)
        grid, edges = dx.histogramdd()
        self.__init__(grid=grid, edges=edges, metadata=self.metadata)

    def _load_plt(self, filename):
        """Initialize Grid from gOpenMol plt file."""
        g = gOpenMol.Plt()
        g.read(filename)
        grid, edges = g.histogramdd()
        self.__init__(grid=grid, edges=edges, metadata=self.metadata)

    def export(self, filename, file_format=None):
        """export density to file using the given format.

        The format can also be deduced from the suffix of the filename
        though the *format* keyword takes precedence.

        The default format for export() is 'dx'.  Use 'dx' for
        visualization.

        Implemented formats:

        dx
            :mod:`OpenDX`
        pickle
            pickle (use :meth:``Grid.load` to restore); :meth:`Grid.save`
            is simpler than ``export(format='python')``.

        """
        exporter = self._get_exporter(filename, file_format=file_format)
        exporter(filename)

    def _export_python(self, filename):
        """Pickle the Grid object

        The object is dumped as a dictionary with grid and edges: This
        is sufficient to recreate the grid object with __init__().
        """
        root, ext = os.path.splitext(filename)
        filename = root + ".pickle"

        data = dict(grid=self.grid, edges=self.edges, metadata=self.metadata)
        with open(filename, 'wb') as f:
            cPickle.dump(data, f, cPickle.HIGHEST_PROTOCOL)

    def _export_dx(self, filename):
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
            'varying, then y, then finally x, i.e. z is the innermost loop.'
        ]

        # write metadata in comments section
        if self.metadata:
            comments.append('Meta data stored with the python Grid object:')
        for k in self.metadata:
            comments.append('   ' + str(k) + ' = ' + str(self.metadata[k]))
        comments.append(
            '(Note: the VMD dx-reader chokes on comments below this line)')

        components = dict(
            positions=OpenDX.gridpositions(1, self.grid.shape, self.origin,
                                           self.delta),
            connections=OpenDX.gridconnections(2, self.grid.shape),
            data=OpenDX.array(3, self.grid), )
        dx = OpenDX.field('density', components=components, comments=comments)
        dx.write(filename)

    def save(self, filename):
        """Save a grid object to <filename>.pickle

           Grid.save(filename)

        Internally, this calls Grid.export(filename,format="python"). A grid can be
        regenerated from the saved data with

           g = Grid(filename=<filename>)

        """
        self.export(filename, file_format="pickle")

    def centers(self):
        """Returns the coordinates of the centers of all grid cells as an
        iterator."""
        for idx in numpy.ndindex(self.grid.shape):
            yield self.delta * numpy.array(idx) + self.origin

    def check_compatible(self, other):
        """Check if *other* can be used in an arithmetic operation.

        1) *other* is a scalar
        2) *other* is a grid defined on the same edges

        :Raises: :exc:`TypeError` if not compatible.
        """
        if not (numpy.isreal(other) or self == other):
            raise TypeError(
                "The argument can not be arithmetically combined with the grid. "
                "It must be a scalar or a grid with identical edges. "
                "Use Grid.resample(other.edges) to make a new grid that is "
                "compatible with other.")
        return True

    def _interpolationFunctionFactory(self, spline_order=None, cval=None):
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
        import scipy.ndimage

        if spline_order is None:
            # must be compatible with whatever :func:`scipy.ndimage.spline_filter` takes.
            spline_order = self.interpolation_spline_order
        if cval is None:
            cval = self.interpolation_cval

        data = self.grid
        if cval is None:
            cval = data.min()
        try:
            # masked arrays, fill with min: should keep spline happy
            _data = data.filled(cval)
        except AttributeError:
            _data = data

        coeffs = scipy.ndimage.spline_filter(_data, order=spline_order)
        x0 = self.origin
        dx = self.delta

        def _transform(cnew, c0, dc):
            return (numpy.atleast_1d(cnew) - c0) / dc

        def interpolatedF(*coordinates):
            """B-spline function over the data grid(x,y,z).

            interpolatedF([x1,x2,...],[y1,y2,...],[z1,z2,...]) -> F[x1,y1,z1],F[x2,y2,z2],...

            Example usage for resampling::
              >>> XX,YY,ZZ = numpy.mgrid[40:75:0.5, 96:150:0.5, 20:50:0.5]
              >>> FF = _interpolationFunction(XX,YY,ZZ)
            """
            _coordinates = numpy.array(
                [_transform(coordinates[i], x0[i], dx[i]) for i in range(len(
                    coordinates))])
            return scipy.ndimage.map_coordinates(coeffs,
                                                 _coordinates,
                                                 prefilter=False,
                                                 mode='nearest',
                                                 cval=cval)
        # mode='wrap' would be ideal but is broken: https://github.com/scipy/scipy/issues/1323
        return interpolatedF

    def __eq__(self, other):
        if not isinstance(other, Grid):
            return False
        return numpy.all(other.grid == self.grid) and numpy.all(
            other.origin == self.origin) and other.edges == self.edges

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        self.check_compatible(other)
        return Grid(self.grid + _grid(other), edges=self.edges)

    def __sub__(self, other):
        self.check_compatible(other)
        return Grid(self.grid - _grid(other), edges=self.edges)

    def __mul__(self, other):
        self.check_compatible(other)
        return Grid(self.grid * _grid(other), edges=self.edges)

    def __div__(self, other):
        self.check_compatible(other)
        return Grid(self.grid / _grid(other), edges=self.edges)

    def __truediv__(self, other):
        self.check_compatible(other)
        return Grid(self.grid / _grid(other), edges=self.edges)

    def __pow__(self, other):
        self.check_compatible(other)
        return Grid(numpy.power(self.grid, _grid(other)), edges=self.edges)

    def __radd__(self, other):
        self.check_compatible(other)
        return Grid(_grid(other) + self.grid, edges=self.edges)

    def __rsub__(self, other):
        self.check_compatible(other)
        return Grid(_grid(other) - self.grid, edges=self.edges)

    def __rmul__(self, other):
        self.check_compatible(other)
        return Grid(_grid(other) * self.grid, edges=self.edges)

    def __rdiv__(self, other):
        self.check_compatible(other)
        return Grid(_grid(other) / self.grid, edges=self.edges)

    def __rtruediv__(self, other):
        return Grid(_grid(other) / self.grid, edges=self.edges)

    def __rpow__(self, other):
        self.check_compatible(other)
        return Grid(numpy.power(_grid(other), self.grid), edges=self.edges)

    def __repr__(self):
        try:
            bins = self.grid.shape
        except AttributeError:
            bins = "no"
        return '<Grid with ' + str(bins) + ' bins>'


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
    lens = list(map(len, arrs))
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz *= s

    ans = []
    for i, arr in enumerate(arrs):
        slc = [1] * dim
        slc[i] = lens[i]
        arr2 = numpy.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j != i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)

    return tuple(ans)
