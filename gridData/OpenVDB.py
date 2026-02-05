r"""
:mod:`~gridData.OpenVDB` --- routines to write OpenVDB files
=============================================================

The `OpenVDB format`_ is used by Blender_ and other VFX software for
volumetric data.

.. _`OpenVDB format`: https://www.openvdb.org
.. _Blender: https://www.blender.org/

This module uses the openvdb_ library to write OpenVDB files.
 
.. _openvdb: https://github.com/AcademySoftwareFoundation/openvdb

.. Note:: This module implements a simple writer for 3D regular grids,
          sufficient to export density data for visualization in Blender_.
          See the `Blender volume docs`_ for details on importing VDB files.
          
.. _`Blender volume docs`: https://docs.blender.org/manual/en/latest/modeling/volumes/introduction.html

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
(File -> Import -> OpenVDB) or (shift+A -> Volume -> Import OpenVDB). See `importing VDB in Blender`_ for details.

.. _`importing VDB in Blender`: https://docs.blender.org/manual/en/latest/modeling/geometry_nodes/input/import/vdb.html


Building an OpenVDB field from a numpy array
---------------------------------------------

If you want to create VDB files without using the Grid class,
you can directly use the OpenVDB field API. This is useful
for custom workflows or when integrating with other libraries.

Requires:

grid
    numpy 3D array
origin
    cartesian coordinates of the center of the (0,0,0) grid cell
delta
    n x n array with the length of a grid cell along each axis

Example::

  from gridData import OpenVDB
  vdb_field = OpenVDB.field('density')
  vdb_field.populate(grid, origin, delta)
  vdb_field.write('output.vdb')


Classes and functions
---------------------

"""

import numpy

try:
    import openvdb as vdb

except ImportError:
    vdb = None


class OpenVDBField(object):
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

    def __init__(self, grid, origin, delta, name="density", tolerance=1e-10):
        """Initialize an OpenVDB field.

        Parameters
        ----------
        grid : numpy.ndarray
            3D numpy array with the data
        origin : numpy.ndarray
            Coordinates of the center of grid cell [0,0,0]
        delta : numpy.ndarray
            Grid spacing (can be 1D array or diagonal matrix)
        name : str
            Name of the grid (will be visible in Blender), default 'density'
        tolerance : float
            Values below this tolerance are treated as background (sparse),
            default 1e-10

        Raises
        ------
        ImportError
            If openvdb is not installed
        ValueError
            If grid is not 3D, or if delta is not 1D/2D or describes
            non-orthorhombic cell

        """
        if vdb is None:
            raise ImportError(
                "openvdb is required to write VDB files. "
                "Install it with: conda install -c conda-forge openvdb"
            )
        self.name = name
        self.tolerance = tolerance
        self._populate(grid, origin, delta)

    def _populate(self, grid, origin, delta):
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
            If grid is not 3D, or if delta is not 1D/2D or describes
            non-orthorhombic cell

        """
        grid = numpy.asarray(grid)
        if grid.ndim != 3:
            raise ValueError(f"OpenVDB only supports 3D grids, got {grid.ndim}D")

        self.grid = grid.astype(numpy.float32)
        self.grid = numpy.ascontiguousarray(self.grid, dtype=numpy.float32)

        self.origin = numpy.asarray(origin)

        # Handle delta: could be 1D array or diagonal matrix
        delta = numpy.asarray(delta)
        if delta.ndim == 2:
            if delta.shape != (3, 3):
                raise ValueError("delta as a matrix must be 3x3")

            if not numpy.allclose(delta, numpy.diag(numpy.diag(delta))):
                raise ValueError("Non-orthorhombic cells are not supported")

            self.delta = numpy.diag(delta)

        elif delta.ndim == 1:
            if len(delta) != 3:
                raise ValueError("delta must have length-3 for 3D grids")
            self.delta = delta

        else:
            raise ValueError(
                "delta must be either a length-3 vector or a 3x3 diagonal matrix"
            )

    def write(self, filename):
        """Write the field to an OpenVDB file.

        Parameters
        ----------
        filename : str
            Output filename (should end in .vdb)

        """

        vdb_grid = vdb.FloatGrid()
        vdb_grid.name = self.name

        vdb_grid.transform = vdb.createLinearTransform()
        vdb_grid.transform.preScale(self.delta.tolist())
        vdb_grid.transform.postTranslate(self.origin.tolist())

        vdb_grid.copyFromArray(self.grid, tolerance=self.tolerance)
        vdb_grid.prune()

        vdb.write(filename, grids=[vdb_grid])
