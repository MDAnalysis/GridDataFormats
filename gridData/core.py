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
# Having consistent truedivision in this module is essential so that
# its behavior is fully consistent in Python 2 and Python 3.
from __future__ import absolute_import, division

import six
from six.moves import cPickle, range, zip

import os
import errno

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
                 metadata={}, interpolation_spline_order=3,
                 file_format=None):
        """
        Create a Grid object from data.

        From a numpy.histogramdd()::

          grid,edges = numpy.histogramdd(...)
          g = Grid(grid,edges=edges)

        From an arbitrary grid::

          g = Grid(grid,origin=origin,delta=delta)

        From a saved file::

          g = Grid(filename)

        or ::

          g = Grid()
          g.load(filename)

        :Arguments:
          grid
            histogram or density, defined on numpy nD array, or filename
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
          file_format
             file format; only necessary when ``grid`` is a filename (see :meth:`Grid.load`);
             default is ``None`` and the file format is autodetected.


        .. versionchanged:: 0.5.0
           New *file_format* keyword argument.

        """
        # file formats are guess from extension == lower case key
        self._exporters = {
            'DX': self._export_dx,
            'PKL': self._export_python,
            'PICKLE': self._export_python,  # compatibility
            'PYTHON': self._export_python,  # compatibility
        }
        self._loaders = {
            'CCP4': self._load_cpp4,
            'DX': self._load_dx,
            'PLT': self._load_plt,
            'PKL': self._load_python,
            'PICKLE': self._load_python,  # compatibility
            'PYTHON': self._load_python,  # compatibility
        }

        self.metadata = metadata
        self.__interpolated = None  # cache for interpolated grid
        self.__interpolation_spline_order = interpolation_spline_order
        self.interpolation_cval = None  # default to using min(grid)

        if grid is not None:
            if isinstance(grid, six.string_types):
                # can probably safely try to load() it...
                filename = grid
            else:
                try:
                    # Can we read this as a file?
                    # Use str(x) to work with py.path.LocalPath and pathlib.Path instances
                    # even for Python < 3.6
                    with open(str(grid), 'rb'):
                        pass
                except (OSError, IOError):
                    # no, this is probably an array-like thingy
                    filename = None
                else:
                    # yes, let's use it as a file
                    filename = str(grid)

            if filename is not None:
                self.load(filename, file_format=file_format)
            else:
                if edges is not None:
                    # set up from histogramdd-type data
                    self.grid = numpy.asanyarray(grid)
                    self.edges = edges
                    self._update()
                elif origin is not None and delta is not None:
                    # setup from generic data
                    origin = numpy.asanyarray(origin)
                    delta = numpy.asanyarray(delta)
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
                    self.grid = numpy.asanyarray(grid)
                    self._update()
                else:
                    raise ValueError("Wrong/missing data to set up Grid. Use Grid() or "
                                     "Grid(grid=<array>, edges=<list>) or "
                                     "Grid(grid=<array>, origin=(x0, y0, z0), delta=(dx, dy, dz)):\n"
                                     "grid={0} edges={1} origin={2} delta={3}".format(
                                         grid, edges, origin, delta))

    @property
    def interpolation_spline_order(self):
        """Order of the B-spline interpolation of the data.

        3 = cubic; 4 & 5 are also supported

        Only choose values that are acceptable to
        :func:`scipy.ndimage.spline_filter`!

        See Also
        --------
        interpolated
        """

        return self.__interpolation_spline_order

    @interpolation_spline_order.setter
    def interpolation_spline_order(self, x):
        """Setting the ``interpolation_spline_order`` updates :func:`interpolated`

        Because we cache the interpolation function, we need to rebuild the
        cache whenever the interpolation order changes: this is
        handled by :meth:`_update`

        """
        self.__interpolation_spline_order = x
        self._update()


    def resample(self, edges):
        """Resample data to a new grid with edges *edges*.

        This method creates a new grid with the data from the current
        grid resampled to a regular grid specified by *edges*.  The
        order of the interpolation is set by
        :attr:`Grid.interpolation_spline_order`: change the value
        *before* calling :meth:`resample`.

        Parameters
        ----------
        edges : tuple of arrays or Grid
             edges of the new grid or a :class:`Grid` instance that
             provides :attr:`Grid.edges`

        Returns
        -------
        Grid
             a new :class:`Grid` with the data interpolated over the
             new grid cells


        Examples
        --------

        Providing *edges* (a tuple of three arrays, indicating the
        boundaries of each grid cell)::

          g = grid.resample(edges)

        As a convenience, one can also supply another :class:`Grid` as
        the argument for this method ::

          g = grid.resample(othergrid)

        and the edges are taken from :attr:`Grid.edges`.

        """
        try:
            edges = edges.edges  # can also supply another Grid
        except AttributeError:
            pass
        midpoints = self._midpoints(edges)
        coordinates = ndmeshgrid(*midpoints)
        # feed a meshgrid to generate all points
        newgrid = self.interpolated(*coordinates)
        return self.__class__(newgrid, edges)

    def resample_factor(self, factor):
        """Resample to a new regular grid.


        Parameters
        ----------
        factor : float
            The number of grid cells are scaled with `factor` in each
            dimension, i.e., ``factor * N_i`` cells along each
            dimension i.


        Returns
        -------
        Grid


        See Also
        --------
        resample

        """
        # new number of edges N' = (N-1)*f + 1
        newlengths = [(N - 1) * float(factor) + 1 for N in self._len_edges()]
        edges = [numpy.linspace(start, stop, num=int(N), endpoint=True)
                 for (start, stop, N) in
                 zip(self._min_edges(), self._max_edges(), newlengths)]
        return self.resample(edges)

    def _update(self):
        """compute/update all derived data

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

        The :func:`interpolated` function allows one to obtain data
        values for any values of the coordinates::

           interpolated([x1,x2,...],[y1,y2,...],[z1,z2,...]) -> F[x1,y1,z1],F[x2,y2,z2],...

        The interpolation order is set in
        :attr:`Grid.interpolation_spline_order`.

        The interpolated function is computed once and is cached for better
        performance. Whenever :attr:`~Grid.interpolation_spline_order` is
        modified, :meth:`Grid.interpolated` is recomputed.

        The value for unknown data is set in :attr:`Grid.interpolation_cval`
        (TODO: also recompute when ``interpolation_cval`` value is changed.)

        Example
        -------
        Example usage for resampling::

           XX, YY, ZZ = numpy.mgrid[40:75:0.5, 96:150:0.5, 20:50:0.5]
           FF = interpolated(XX, YY, ZZ)

        Note
        ----
        Values are interpolated with a spline function. It is possible
        that the spline will generate values that would not normally
        appear in the data. For example, a density is non-negative but
        a cubic spline interpolation can generate negative values,
        especially at the boundary between 0 and high values.

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
        filename = str(filename)
        if not os.path.exists(filename):
            # check before we try to detect the file type because
            # _guess_fileformat() does not work well with things that
            # are not really a file
            raise IOError(errno.ENOENT, "file not found", filename)
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

    def export(self, filename, file_format=None, type=None, typequote='"'):
        """export density to file using the given format.

        The format can also be deduced from the suffix of the filename
        though the *format* keyword takes precedence.

        The default format for export() is 'dx'.  Use 'dx' for
        visualization.

        Implemented formats:

        dx
            :mod:`OpenDX`
        pickle
            pickle (use :meth:`Grid.load` to restore); :meth:`Grid.save`
            is simpler than ``export(format='python')``.

        Parameters
        ----------

        filename : str
            name of the output file

        file_format : {'dx', 'pickle', None} (optional)
            output file format, the default is "dx"

        type : str (optional)
            for DX, set the output DX array type, e.g., "double" or "float".
            By default (``None``), the DX type is determined from the numpy
            dtype of the array of the grid (and this will typically result in
            "double").

            .. versionadded:: 0.4.0

        typequote : str (optional)
            For DX, set the character used to quote the type string;
            by default this is a double-quote character, '"'.
            Custom parsers like the one from NAMD-GridForces (backend for MDFF)
            expect no quotes, and typequote='' may be used to appease them.

            .. versionadded:: 0.5.0

        """
        filename = str(filename)
        exporter = self._get_exporter(filename, file_format=file_format)
        exporter(filename, type=type, typequote=typequote)

    # note: the _export_FORMAT() methods all take the filename as a mandatory
    # argument. They can process kwargs but they are not required to do
    # so. However, they must ignore any kwargs that they are not processing.

    def _export_python(self, filename, **kwargs):
        """Pickle the Grid object

        The object is dumped as a dictionary with grid and edges: This
        is sufficient to recreate the grid object with __init__().
        """
        data = dict(grid=self.grid, edges=self.edges, metadata=self.metadata)
        with open(filename, 'wb') as f:
            cPickle.dump(data, f, cPickle.HIGHEST_PROTOCOL)

    def _export_dx(self, filename, type=None, typequote='"', **kwargs):
        """Export the density grid to an OpenDX file.

        The file format is the simplest regular grid array and it is
        also understood by VMD's and Chimera's DX reader; PyMOL
        requires the dx `type` to be set to "double".

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
            data=OpenDX.array(3, self.grid, type=type, typequote=typequote),
        )
        dx = OpenDX.field('density', components=components, comments=comments)
        dx.write(filename)

    def save(self, filename):
        """Save a grid object to <filename>.pickle

        Internally, this calls
        ``Grid.export(filename, format="python")``. A grid can be
        regenerated from the saved data with ::

           g = Grid(filename="grid.pickle")

        .. note::
           The pickle format depends on the Python version and
           therefore it is not guaranteed that a grid saved with, say,
           Python 2.7 can also be read with Python 3.5. The OpenDX format
           is a better alternative for portability.

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
        return numpy.all(other.grid == self.grid) and \
            numpy.all(other.origin == self.origin) and \
            numpy.all(numpy.all(other_edge == self_edge) for other_edge, self_edge in
                      zip(other.edges, self.edges))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid + _grid(other), edges=self.edges)

    def __sub__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid - _grid(other), edges=self.edges)

    def __mul__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid * _grid(other), edges=self.edges)

    def __truediv__(self, other):
        # truediv will always do true division (in Python 2 and Python 3);
        # we use from __future__ include division everywhere
        self.check_compatible(other)
        return self.__class__(self.grid / _grid(other), edges=self.edges)

    def __div__(self, other):
        # in Python 2 only (without __future__.division): will do "classic division"
        # https://docs.python.org/2/reference/datamodel.html#object.__div__
        if not six.PY2:
            raise NotImplementedError("__div__ is only available in Python 2, use __truediv__")
        self.check_compatible(other)
        return self.__class__(self.grid.__div__(_grid(other)), edges=self.edges)

    def __floordiv__(self, other):
        self.check_compatible(other)
        return self.__class__(self.grid // _grid(other), edges=self.edges)

    def __pow__(self, other):
        self.check_compatible(other)
        return self.__class__(numpy.power(self.grid, _grid(other)), edges=self.edges)

    def __radd__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) + self.grid, edges=self.edges)

    def __rsub__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) - self.grid, edges=self.edges)

    def __rmul__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) * self.grid, edges=self.edges)

    def __rtruediv__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) / self.grid, edges=self.edges)

    def __rdiv__(self, other):
        # in Python 2 only (without __future__.division): will do "classic division"
        # https://docs.python.org/2/reference/datamodel.html#object.__div__
        if not six.PY2:
            raise NotImplementedError("__rdiv__ is only available in Python 2, use __rtruediv__")
        self.check_compatible(other)
        return self.__class__(self.grid.__rdiv__(_grid(other)), edges=self.edges)

    def __rfloordiv__(self, other):
        self.check_compatible(other)
        return self.__class__(_grid(other) // self.grid, edges=self.edges)

    def __rpow__(self, other):
        self.check_compatible(other)
        return self.__class__(numpy.power(_grid(other), self.grid), edges=self.edges)

    def __repr__(self):
        try:
            bins = self.grid.shape
        except AttributeError:
            bins = "no"
        return '<{0} with {1!r} bins>'.format(self.__class__, bins)


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
        arr2 = numpy.asanyarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j != i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)

    return tuple(ans)
