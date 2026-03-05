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

  from gridData import Grid
  
  g = Grid("data.dx")
  g.export("data.vdb")


Classes and functions
---------------------

"""

import numpy
import warnings
from dataclasses import dataclass

try:
    import openvdb as vdb

except ImportError:
    vdb = None


@dataclass
class DownCastTo:
    gridType: str


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

      import gridData.OpenVDB as OpenVDB

      vdb_field = OpenVDB.OpenVDBField('density')
      vdb_field.write('output.vdb')

    Or use directly from Grid::

      g = Grid(...)
      g.export('output.vdb', format='vdb')

    """

    def __init__(
        self,
        grid=None,
        origin=None,
        delta=None,
        name="density",
        tolerance=None,
        metadata=None,
    ):
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
        tolerance : float (optional)
            Values below this tolerance are treated as background (sparse),
            default None
        metadata : dict (optional)
            Additional metadata to embed in the VDB file.

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
        if metadata is not None:
            self.metadata = metadata
        else:
            self.metadata = {}

        if grid is not None:
            if isinstance(grid, vdb.GridBase):
                self.vdb_grid = grid
                print("hello")
                self._extract_from_vdb_grid()
            else:
                self._populate(grid, origin, delta)
                self.vdb_grid = self._create_openvdb_grid()
        else:
            self.grid = None
            self.origin = None
            self.delta = None
            self.vdb_grid = None

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
        self.grid = numpy.asarray(grid)
        if self.grid.ndim != 3:
            raise ValueError(f"OpenVDB only supports 3D grids, got {grid.ndim}D")

        self.grid = numpy.ascontiguousarray(self.grid)

        if origin is None:
            raise ValueError("origin must be provided for VDB export")

        if len(origin) != 3:
            raise ValueError("origin must be length-3")

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

    def _get_best_grid_type(self):
        """Selects the suitable OpenVDB grid type

        Returns
        -------
        openvdb.GridBase

        Raises
        ------
        TypeError
            If dtype is not supported or no suitable grid type is available
        """
        datatypes = {
            numpy.dtype("bool"): ["BoolGrid"],
            numpy.dtype("int8"): ["Int32Grid", "FloatGrid"],
            numpy.dtype("uint8"): ["Int32Grid", "FloatGrid"],
            numpy.dtype("int16"): ["Int32Grid", "FloatGrid"],
            numpy.dtype("uint16"): ["Int32Grid", "FloatGrid"],
            numpy.dtype("int32"): ["Int32Grid", DownCastTo("FloatGrid")],
            numpy.dtype("uint32"): [DownCastTo("Int32Grid"), DownCastTo("FloatGrid")],
            numpy.dtype("int64"): ["Int64Grid", DownCastTo("FloatGrid")],
            numpy.dtype("uint64"): ["Int64Grid", DownCastTo("FloatGrid")],
            numpy.dtype("float16"): ["HalfGrid", "FloatGrid"],
            numpy.dtype("float32"): ["FloatGrid"],
            numpy.dtype("float64"): ["DoubleGrid", DownCastTo("FloatGrid")],
        }

        try:
            vdb_gridtypes = datatypes[self.grid.dtype]
        except KeyError:
            raise TypeError(f"Data type {self.grid.dtype} not supported for VDB")

        VDB_Grid = None
        selected_gridtype = None
        is_downcast = False

        for gridtype in vdb_gridtypes:
            if isinstance(gridtype, DownCastTo):
                gridtype_name = gridtype.gridType
                is_downcast = True
            else:
                gridtype_name = gridtype
                is_downcast = False

            try:
                VDB_Grid = getattr(vdb, gridtype_name)
            except AttributeError:
                continue
            else:
                selected_gridtype = gridtype_name
                break
        else:
            raise TypeError(
                f"Could not find any VDB grid type for numpy dtype {self.grid.dtype}"
            )

        if is_downcast:
            warnings.warn(
                f"Grid type {vdb_gridtypes[0]} not available. Using {selected_gridtype} instead. Data may lose precision.",
                RuntimeWarning,
            )
        return VDB_Grid()

    def _create_openvdb_grid(self):
        """Create and populate an OpenVDB grid

        Returns
        -------
        openvdb.GridBase

        """

        vdb_grid = self._get_best_grid_type()

        vdb_grid.name = self.name

        vdb_grid.transform = vdb.createLinearTransform()
        vdb_grid.transform.preScale(self.delta.tolist())
        vdb_grid.transform.postTranslate(self.origin.tolist())

        if self.metadata:
            for key, val in self.metadata.items():
                try:
                    vdb_grid[key] = val
                except (TypeError, ValueError) as e:
                    warnings.warn(f"Could not set metadata '{key}': {e}", UserWarning)

        if isinstance(vdb_grid, vdb.BoolGrid) and (
            self.tolerance is None or self.tolerance == 0
        ):
            vdb_grid.copyFromArray(self.grid)
            vdb_grid.prune(tolerance=False)
        else:
            if self.tolerance is None:
                vdb_grid.copyFromArray(self.grid)
            else:
                vdb_grid.copyFromArray(self.grid, tolerance=self.tolerance)

            vdb_grid.prune()

        return vdb_grid

    def _extract_from_vdb_grid(self):
        """Extract numpy array, origin, delta from stored VDB grid.

        This method converts the sparse VDB grid to a dense numpy array
        and extracts the transform information.
        """
        if self.vdb_grid.name:
            self.name = self.vdb_grid.name

        for key in self.vdb_grid.metadata:
            try:
                self.metadata[key] = self.vdb_grid[key]
            except (TypeError, ValueError):
                pass

        transformation = self.vdb_grid.transform

        v_origin = numpy.array(transformation.indexToWorld([0, 0, 0]))
        v_delta = numpy.array(transformation.indexToWorld([1, 1, 1])) - v_origin

        self.origin = v_origin
        self.delta = v_delta

        vdb_native_dtype = {
            "BoolGrid": numpy.dtype("bool"),
            "Int32Grid": numpy.dtype("int32"),
            "Int64Grid": numpy.dtype("int64"),
            "FloatGrid": numpy.dtype("float32"),
            "DoubleGrid": numpy.dtype("float64"),
            "HalfGrid": numpy.dtype("float16"),
            "Vec3SGrid": numpy.dtype("float32"),
            "Vec3DGrid": numpy.dtype("float64"),
        }
        dtype = vdb_native_dtype.get(
            type(self.vdb_grid).__name__, numpy.dtype("float32")
        )

        bbox = self.vdb_grid.evalActiveVoxelBoundingBox()

        if bbox is None or bbox[0] == bbox[1]:
            self.grid = numpy.zeros((0, 0, 0), dtype=dtype)
            return

        shape = tuple(numpy.array(bbox[1]) - numpy.array(bbox[0]) + 1)

        self.grid = numpy.zeros(shape, dtype=dtype)
        print(dtype)
        self.vdb_grid.copyToArray(self.grid, ijk=bbox[0])

        if not numpy.all(numpy.array(bbox[0]) == 0):
            self.origin = numpy.array(
                transformation.indexToWorld(numpy.array(bbox[0]).tolist())
            )

    def write(self, filename):
        """Write the field to an OpenVDB file.

        Parameters
        ----------
        filename : str
            Output filename (should end in .vdb)
        """
        if self.vdb_grid is None:
            raise ValueError("No grid data to write")
        vdb.write(filename, grids=[self.vdb_grid])
