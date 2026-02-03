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

        outfile = str(tmpdir / "metadata.vdb")
        g.export(outfile)

        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == "test_density"

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

    def test_write_vdb_zero_threshold(self, tmpdir):
        data = np.ones((3, 3, 3), dtype=np.float32) * 1e-11
        data[1, 1, 1] = 1.0

        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "threshold.vdb")
        g.export(outfile)
        assert tmpdir.join("threshold.vdb").exists()

        grids, metadata = vdb.readAll(outfile)
        grid_vdb = grids[0]
        acc = grid_vdb.getAccessor()

        assert acc.getValue((1, 1, 1)) == pytest.approx(data[1, 1, 1])

        val, is_active = grid_vdb.getConstAccessor().probeValue((0, 0, 0))

        assert (not is_active) or (val == pytest.approx(0.0))

    def test_vdb_non_orthrhombic_raises(self):
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


def test_vdb_import_error():
    with patch("gridData.OpenVDB.vdb", None):
        with pytest.raises(ImportError, match="openvdb is required"):
            gridData.OpenVDB.OpenVDBField(
                np.ones((3, 3, 3)), origin=[0, 0, 0], delta=[1, 1, 1]
            )
