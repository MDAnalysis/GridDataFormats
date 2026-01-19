r"""
:mod:`~gridData.OpenVDB` --- routines to write OpenVDB files
=============================================================

The OpenVDB format is used by Blender and other VFX software for
volumetric data. See https://www.openvdb.org

pyopenvdb: https://github.com/AcademySoftwareFoundation/openvdb

.. Note:: This module implements a simple writer for 3D regular grids,
          sufficient to export density data for visualization in Blender.

The OpenVDB format uses a sparse tree structure to efficiently store
volumetric data. It is the native format for Blender's volume system.


Writing OpenVDB files
---------------------

If you have a :class:`~gridData.core.Grid` object, you can write it to
OpenVDB format::

  from gridData import Grid
  g = Grid("data.dx")
  g.export("data.vdb")

This will create a file that can be imported directly into Blender
(File -> Import -> OpenVDB).


Building an OpenVDB field from a numpy array
---------------------------------------------

Requires:

grid
    numpy 3D array
origin
    cartesian coordinates of the center of the (0,0,0) grid cell
delta
    n x n array with the length of a grid cell along each axis

Example::

  import OpenVDB
  vdb_field = OpenVDB.field('density')
  vdb_field.populate(grid, origin, delta)
  vdb_field.write('output.vdb')


Classes and functions
---------------------

"""

import numpy
import warnings

try:
    import pyopenvdb as vdb
except ImportError:
    try:
        import openvdb as vdb
    except ImportError:
        vdb = None


class field(object):
    """OpenVDB field object for writing volumetric data.

    This class provides a simple interface to write 3D grid data to
    OpenVDB format, which can be imported into Blender and other
    VFX software.

    The field object holds grid data and metadata, and can write it
    to a .vdb file.

    Example
    -------
    Create a field and write it::

      vdb_field = OpenVDB.field('density')
      vdb_field.populate(grid, origin, delta)
      vdb_field.write('output.vdb')

    Or use directly from Grid::

      g = Grid(...)
      g.export('output.vdb', format='vdb')

    """

    def __init__(self, name='density'):
        """Initialize an OpenVDB field.

        Parameters
        ----------
        name : str
            Name of the grid (will be visible in Blender)

        """
        if vdb is None:
            raise ImportError(
                "pyopenvdb is required to write VDB files. "
                "Install it with: conda install -c conda-forge openvdb"
            )
        self.name = name
        self.grid = None
        self.origin = None
        self.delta = None

    def populate(self, grid, origin, delta):
        """Populate the field with grid data.

        Parameters
        ----------
        grid : numpy.ndarray
            3D numpy array with the data
        origin : numpy.ndarray
            Coordinates of the center of grid cell [0,0,0]
        delta : numpy.ndarray
            Grid spacing (can be 1D array or diagonal matrix)

        Raises
        ------
        ValueError
            If grid is not 3D

        """
        grid = numpy.asarray(grid)
        if grid.ndim != 3:
            raise ValueError(
                "OpenVDB only supports 3D grids, got {}D".format(grid.ndim))

        self.grid = numpy.transpose(grid, (2, 1, 0)).astype(numpy.float32)
        self.origin = numpy.asarray(origin)

        # Handle delta: could be 1D array or diagonal matrix
        delta = numpy.asarray(delta)
        if delta.ndim == 2:
            self.delta = numpy.array([delta[i, i] for i in range(3)])
        else:
            self.delta = delta

    def write(self, filename):
        """Write the field to an OpenVDB file.

        Parameters
        ----------
        filename : str
            Output filename (should end in .vdb)

        """
        if self.grid is None:
            raise ValueError("No data to write. Use populate() first.")

        # Create OpenVDB grid
        vdb_grid = vdb.FloatGrid()
        vdb_grid.name = self.name

        # this is an explicit linear transform using per-axis voxel sizes
        # world = diag(delta) * index + corner_origin
        corner_origin = (self.origin - 0.5 * self.delta).astype(float)

        # Constructing 4x4 row-major matrix where the last row is the translation
        matrix = [
            [float(self.delta[0]), 0.0, 0.0, 0.0],
            [0.0, float(self.delta[1]), 0.0, 0.0],
            [0.0, 0.0, float(self.delta[2]), 0.0],
            [float(corner_origin[0]), float(corner_origin[1]), float(corner_origin[2]), 1.0]
        ]

        transform = vdb.createLinearTransform(matrix)
        vdb_grid.transform = transform

        vdb_grid.background = 0.0

        # Populate the grid
        accessor = vdb_grid.getAccessor()
        threshold = 1e-10 

        mask = numpy.abs(self.grid) > threshold
        indices = numpy.argwhere(mask)
        
        for idx in indices:
            i, j, k = idx
            value = float(self.grid[i, j, k])
            accessor.setValueOn((i, j, k), value)

        vdb.write(filename, grids=[vdb_grid])