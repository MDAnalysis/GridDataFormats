import numpy as np
from numpy.testing import assert_allclose

import pytest
from unittest.mock import patch

import gridData.OpenVDB
from gridData import Grid

from . import datafiles

try:
    import openvdb as vdb

    HAS_OPENVDB = True
except ImportError:
    HAS_OPENVDB = False


@pytest.fixture
def grid345():
    data = np.arange(1, 61, dtype=np.float32).reshape((3, 4, 5))
    g = Grid(data.copy(), origin=np.zeros(3), delta=np.ones(3))
    return data, g


@pytest.mark.skipif(not HAS_OPENVDB, reason="openvdb not installed")
class TestVDBWrite:
    def test_write_vdb_from_grid(self, tmpdir, grid345):
        data, g = grid345

        outfile = str(tmpdir / "test.vdb")
        g.export(outfile, file_format="VDB")

        assert tmpdir.join("test.vdb").exists()

        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1
        assert grids[0].name == "density"

        grid_vdb = grids[0]
        acc = grid_vdb.getAccessor()

        assert grid_vdb.evalActiveVoxelDim() == data.shape

        corners = [
            (0, 0, 0),
            (data.shape[0] - 1, data.shape[1] - 1, data.shape[2] - 1),
            (1, 2, 3),
        ]
        for i, j, k in corners:
            got = acc.getValue((i, j, k))
            assert got == pytest.approx(float(data[i, j, k]))

    def test_write_vdb_default_grid_name(self, tmpdir, grid345):
        data, g = grid345
        g.metadata = {}

        outfile = str(tmpdir / "default_name.vdb")
        g.export(outfile)
        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == "density"

    def test_write_vdb_autodetect_extension(self, tmpdir, grid345):
        data, g = grid345

        outfile = str(tmpdir / "auto.vdb")
        g.export(outfile)

        assert tmpdir.join("auto.vdb").exists()

    def test_write_vdb_with_metadata(self, tmpdir, grid345):
        data, g = grid345
        g.metadata["name"] = "test_density"
        g.metadata["source"] = "test"
        g.metadata["factor"] = 60

        outfile = str(tmpdir / "metadata.vdb")
        g.export(outfile)

        metadata = vdb.readAllGridMetadata(outfile)
        assert metadata[0]["name"] == "test_density"
        assert metadata[0]["source"] == "test"
        assert metadata[0]["factor"] == 60

    def test_write_vdb_origin_and_spacing(self, tmpdir):
        data = np.ones((4, 4, 4), dtype=np.float32)
        origin = np.array([10.0, 20.0, 30.0])
        delta = np.array([0.5, 0.5, 0.5])

        g = Grid(data, origin=origin, delta=delta)
        outfile = str(tmpdir / "transform.vdb")
        g.export(outfile)

        grids, metadata = vdb.readAll(outfile)
        grid_vdb = grids[0]

        voxel_size = grid_vdb.transform.voxelSize()

        spacing = [voxel_size[0], voxel_size[1], voxel_size[2]]
        assert_allclose(spacing, delta, rtol=1e-5)

    def test_write_vdb_from_ccp4(self, tmpdir):
        g = Grid(datafiles.CCP4)
        outfile = str(tmpdir / "from_ccp4.vdb")

        g.export(outfile, file_format="VDB")

        assert tmpdir.join("from_ccp4.vdb").exists()
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1

    def test_vdb_field_direct(self, tmpdir):
        data = np.arange(27).reshape(3, 3, 3).astype(np.float32)

        vdb_field = gridData.OpenVDB.OpenVDBField(
            data, origin=[0, 0, 0], delta=[1, 1, 1], name="direct_test"
        )

        outfile = str(tmpdir / "direct.vdb")
        vdb_field.write(outfile)

        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == "direct_test"

        grid_vdb = grids[0]
        assert grid_vdb.evalActiveVoxelDim() == data.shape
        acc = grid_vdb.getAccessor()

        assert acc.getValue((0, 0, 0)) == pytest.approx(float(data[0, 0, 0]))
        assert acc.getValue((1, 1, 1)) == pytest.approx(float(data[1, 1, 1]))
        assert acc.getValue((2, 2, 2)) == pytest.approx(float(data[2, 2, 2]))

    def test_vdb_field_2d_raises(self):
        data_2d = np.arange(12).reshape(3, 4)

        with pytest.raises(ValueError, match="3D grids"):
            gridData.OpenVDB.OpenVDBField(data_2d, origin=[0, 0], delta=[1, 1])

    def test_write_vdb_nonuniform_spacing(self, tmpdir):
        data = np.ones((3, 3, 3), dtype=np.float32)
        delta = np.array([0.5, 1.0, 1.5])
        g = Grid(data, origin=[0, 0, 0], delta=delta)

        outfile = str(tmpdir / "nonuniform.vdb")
        g.export(outfile)
        assert tmpdir.join("nonuniform.vdb").exists()

        grids, metadata = vdb.readAll(outfile)
        grid_vdb = grids[0]
        voxel_size = grid_vdb.transform.voxelSize()

        spacing = [voxel_size[0], voxel_size[1], voxel_size[2]]

        assert_allclose(spacing, delta, rtol=1e-5)

    def test_write_vdb_with_delta_matrix(self, tmpdir):
        data = np.ones((3, 3, 3), dtype=np.float32)
        delta = np.diag([1.0, 2.0, 3.0])

        vdb_field = gridData.OpenVDB.OpenVDBField(
            data, origin=[0, 0, 0], delta=delta, name="matrix_delta"
        )

        outfile = str(tmpdir / "matrix_delta.vdb")
        vdb_field.write(outfile)
        assert tmpdir.join("matrix_delta.vdb").exists()

        grids, metadata = vdb.readAll(outfile)
        grid_vdb = grids[0]
        voxel_size = grid_vdb.transform.voxelSize()

        spacing = [voxel_size[0], voxel_size[1], voxel_size[2]]

        assert_allclose(spacing, np.diag(delta), rtol=1e-5)

    def test_write_vdb_sparse_data(self, tmpdir):
        data = np.zeros((10, 10, 10), dtype=np.float32)
        data[2, 3, 4] = 5.0
        data[7, 8, 9] = 10.0

        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "sparse.vdb")
        g.export(outfile)

        assert tmpdir.join("sparse.vdb").exists()
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1

        grid_vdb = grids[0]
        acc = grid_vdb.getAccessor()
        assert acc.getValue((2, 3, 4)) == pytest.approx(data[2, 3, 4])
        assert acc.getValue((7, 8, 9)) == pytest.approx(data[7, 8, 9])

    def test_write_vdb_with_zero_tolerance(self, tmpdir):
        data = np.zeros((3, 3, 3), dtype=np.float32)
        data[0, 0, 0] = 1.0
        data[1, 1, 1] = 2.0
        data[2, 2, 2] = 1e-7

        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "zero_tolerance.vdb")

        g.export(outfile, tolerance=0)

        assert tmpdir.join("zero_tolerance.vdb").exists()

        grids, _ = vdb.readAll(outfile)
        grid_vdb = grids[0]
        acc = grid_vdb.getAccessor()

        assert acc.getValue((0, 0, 0)) == pytest.approx(1.0)
        assert acc.getValue((1, 1, 1)) == pytest.approx(2.0)
        assert acc.getValue((2, 2, 2)) == pytest.approx(1e-7)

        outfile_pruned = str(tmpdir / "with_tolerance.vdb")
        g.export(outfile_pruned, tolerance=1e-6)
        grids_pruned, _ = vdb.readAll(outfile_pruned)
        assert not grids_pruned[0].getAccessor().isValueOn((2, 2, 2))

    def test_vdb_non_orthorhombic_raises(self):
        data = np.ones((3, 3, 3), dtype=np.float32)
        delta = np.array(
            [
                [1, 0.1, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )

        with pytest.raises(ValueError, match="Non-orthorhombic"):
            gridData.OpenVDB.OpenVDBField(data, origin=[0, 0, 0], delta=delta)

    def test_delta_matrix_wrong_shape_raises(self):
        data = np.ones((3, 3, 3), dtype=np.float32)
        origin = [0.0, 0.0, 0.0]
        bad_delta = np.eye(2)

        with pytest.raises(ValueError, match="delta as a matrix must be 3x3"):
            gridData.OpenVDB.OpenVDBField(data, origin, bad_delta)

    def test_delta_scalar_raises(self):
        data = np.ones((3, 3, 3), dtype=np.float32)
        origin = [0.0, 0.0, 0.0]
        bad_delta = np.array(1.0)

        with pytest.raises(
            ValueError,
            match="delta must be either a length-3 vector or a 3x3 diagonal matrix",
        ):
            gridData.OpenVDB.OpenVDBField(data, origin, bad_delta)

    def test_delta_1d_wrong_length_raises(self):
        data = np.ones((3, 3, 3), dtype=np.float32)
        origin = [0.0, 0.0, 0.0]
        bad_delta = np.array([1.0, 2.0])

        with pytest.raises(ValueError, match="must have length-3"):
            gridData.OpenVDB.OpenVDBField(data, origin, bad_delta)

    def test_write_vdb_boolean_grid(self, tmpdir):
        data = np.random.random((10, 10, 10)) > 0.5

        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        g.metadata["name"] = "boolean_mask"

        outfile = str(tmpdir / "boolean.vdb")
        g.export(outfile)
        assert tmpdir.join("boolean.vdb").exists()

        grids, _ = vdb.readAll(outfile)
        assert len(grids) == 1

        grid_vdb = grids[0]

        assert isinstance(grid_vdb, vdb.BoolGrid)

        acc = grid_vdb.getAccessor()
        assert acc.getValue((0, 0, 0)) == bool(data[0, 0, 0])
        assert acc.getValue((5, 5, 5)) == bool(data[5, 5, 5])

    def test_write_vdb_with_tolerance_parameter(self, tmpdir, grid345):
        data, g = grid345
        outfile = str(tmpdir / "with_tolerance.vdb")
        g.export(outfile, tolerance=0.1)

        assert tmpdir.join("with_tolerance.vdb").exists()

        grids, _ = vdb.readAll(outfile)
        grid_vdb = grids[0]

        acc = grid_vdb.getAccessor()
        assert acc.getValue((0, 0, 0)) == pytest.approx(float(data[0, 0, 0]))
        assert acc.getValue((1, 2, 3)) == pytest.approx(float(data[1, 2, 3]))

    def test_write_vdb_float64_conversion_warning(self, tmpdir):
        data = np.random.random((5, 5, 5))

        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "float64.vdb")

        with pytest.warns(
            RuntimeWarning, match="Grid type DoubleGrid not available.*FloatGrid"
        ):
            g.export(outfile)

        grids, _ = vdb.readAll(outfile)
        assert isinstance(grids[0], vdb.FloatGrid)

    def test_write_vdb_int32_conversion_warning(self, tmpdir):
        data = np.arange(27, dtype=np.int32).reshape((3, 3, 3))

        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "int32.vdb")

        with pytest.warns(
            RuntimeWarning, match="Grid type Int32Grid not available.*FloatGrid"
        ):
            g.export(outfile)

        grids, _ = vdb.readAll(outfile)
        assert isinstance(grids[0], vdb.FloatGrid)

        acc = grids[0].getAccessor()
        assert acc.getValue((0, 0, 0)) == pytest.approx(float(data[0, 0, 0]))
        assert acc.getValue((1, 1, 1)) == pytest.approx(float(data[1, 1, 1]))

    def test_write_vdb_unsupported_dtype_raises(self):
        data_complex = np.ones((3, 3, 3), dtype=np.complex64)

        with pytest.raises(TypeError, match="Data type.*not supported for VDB"):
            gridData.OpenVDB.OpenVDBField(
                data_complex, origin=[0, 0, 0], delta=[1, 1, 1]
            )

    def test_openvdb_field_empty_initialization(self, tmpdir):
        vdb_field = gridData.OpenVDB.OpenVDBField()

        assert vdb_field.grid is None
        assert vdb_field.origin is None
        assert vdb_field.delta is None
        assert vdb_field.vdb_grid is None
        assert vdb_field.name == "density"
        assert vdb_field.metadata == {}

        data = np.ones((3, 3, 3), dtype=np.float32)
        vdb_field._populate(data, [0, 0, 0], [1, 1, 1])
        vdb_field.vdb_grid = vdb_field._create_openvdb_grid()

        assert vdb_field.grid is not None
        assert vdb_field.origin is not None
        assert vdb_field.delta is not None
        assert vdb_field.vdb_grid is not None

        outfile = str(tmpdir / "empty_init.vdb")
        vdb_field.write(outfile)

        assert tmpdir.join("empty_init.vdb").exists()

    def test_from_grid_openvdb_field(self, grid345):
        data, g = grid345
        vdb_field = gridData.OpenVDB.OpenVDBField.from_grid(g)

        assert vdb_field.name == "density"
        assert vdb_field.vdb_grid is not None
        assert isinstance(vdb_field.native, vdb.FloatGrid)

    def test_convert_to_vdb_native_grid(self, grid345):
        data, g = grid345
        native = g.convert_to("vdb")

        assert isinstance(native, vdb.GridBase)
        assert native.evalActiveVoxelDim() == data.shape

    def test_from_native_grid_shape_values_and_dimension(self, grid345):
        data, g = grid345
        vdb_field = gridData.OpenVDB.OpenVDBField.from_grid(g)

        native_grid = vdb_field.native

        assert native_grid.evalActiveVoxelDim() == data.shape

        acc = native_grid.getAccessor()
        assert acc.getValue((0, 0, 0)) == pytest.approx(float(data[0, 0, 0]))
        assert acc.getValue((1, 2, 3)) == pytest.approx(float(data[1, 2, 3]))

        voxel_size = native_grid.transform.voxelSize()
        assert_allclose(
            [voxel_size[0], voxel_size[1], voxel_size[2]], g.delta, rtol=1e-5
        )

        world = native_grid.transform.indexToWorld((0, 0, 0))
        assert_allclose([world[0], world[1], world[2]], g.origin, rtol=1e-5)

    def test_file_roundtrip_native_vdbgrid(self, tmpdir, grid345):
        data, g = grid345
        g.metadata["name"] = "new_density"

        outfile = str(tmpdir / "roundtrip.vdb")
        g.export(outfile)

        grids, _ = vdb.readAll(outfile)
        assert len(grids) == 1

        grid_vdb = grids[0]
        new_vdb_grid = Grid(grid=grid_vdb)
        assert_allclose(new_vdb_grid.grid, g.grid, rtol=1e-5)
        assert_allclose(new_vdb_grid.origin, g.origin, rtol=1e-5)
        assert_allclose(new_vdb_grid.delta, g.delta, rtol=1e-5)
        assert new_vdb_grid.metadata["name"] == "new_density"
        assert new_vdb_grid.grid.dtype == np.dtype("float32")
        assert_allclose(new_vdb_grid.grid, data, rtol=1e-5)

    def test_extract_from_vdb_grid(self, grid345):
        data, g = grid345
        g.metadata["name"] = "new_density"

        native = g.convert_to("vdb")
        new_vdb_grid = Grid(grid=native)

        assert_allclose(new_vdb_grid.grid, g.grid, rtol=1e-5)
        assert_allclose(new_vdb_grid.origin, g.origin, rtol=1e-5)
        assert_allclose(new_vdb_grid.delta, g.delta, rtol=1e-5)
        assert new_vdb_grid.metadata["name"] == "new_density"
        assert new_vdb_grid.grid.dtype == np.dtype("float32")
        assert_allclose(new_vdb_grid.grid, data, rtol=1e-5)

    @pytest.mark.parametrize(
        "vdb_type,np_dtype", list(gridData.OpenVDB.OpenVDBField._DTYPES_VDB2NP.items())
    )
    def test_extract_dtype_roundtrip(self, vdb_type, np_dtype):
        vdb_cls = getattr(vdb, vdb_type, None)
        if vdb_cls is None:
            pytest.skip(f"{vdb_type} not available in this openvdb build")

        native = vdb_cls()
        native.name = "test"
        acc = native.getAccessor()
        acc.setValueOn((0, 0, 0), True if np_dtype == np.dtype("bool") else 1)
        acc.setValueOn((1, 1, 1), True if np_dtype == np.dtype("bool") else 2)

        field = gridData.OpenVDB.OpenVDBField(grid=native)

        assert field.grid.dtype == np_dtype
        assert field.grid.shape != (0, 0, 0)

    def test_extract_unknown_vdb_type_warns(self):
        native = vdb.FloatGrid()
        acc = native.getAccessor()
        acc.setValueOn((0, 0, 0), 1.0)

        patched = {
            k: v
            for k, v in gridData.OpenVDB.OpenVDBField._DTYPES_VDB2NP.items()
            if k != "FloatGrid"
        }
        with patch.object(gridData.OpenVDB.OpenVDBField, "_DTYPES_VDB2NP", patched):
            with pytest.warns(RuntimeWarning, match="Unknown VDB grid type"):
                field = gridData.OpenVDB.OpenVDBField(grid=native)

        assert field.grid.dtype == np.dtype("float32")

    def test_extract_empty_vdb_grid_warns(self):
        native = vdb.FloatGrid()
        native.clear()

        with pytest.warns(RuntimeWarning, match="no active voxels"):
            field = gridData.OpenVDB.OpenVDBField(grid=native)

        assert field.grid.shape == (0, 0, 0)


@pytest.mark.skipif(
    not HAS_OPENVDB, reason="Need openvdb to test import error handling"
)
def test_vdb_import_error():
    with patch("gridData.OpenVDB.vdb", None):
        with pytest.raises(ImportError, match="openvdb is required"):
            gridData.OpenVDB.OpenVDBField(
                np.ones((3, 3, 3)), origin=[0, 0, 0], delta=[1, 1, 1]
            )
