import numpy as np
from numpy.testing import assert_allclose, assert_equal

import pytest

import gridData.OpenVDB
from gridData import Grid

from . import datafiles


try:
    import pyopenvdb as vdb
    HAS_OPENVDB = True
except ImportError:
    try:
        import openvdb as vdb
        HAS_OPENVDB = True
    except ImportError:
        HAS_OPENVDB = False


@pytest.mark.skipif(not HAS_OPENVDB, reason="pyopenvdb/openvdb not installed")
class TestVDBWrite:
    """Test OpenVDB file format writing"""

    def test_write_vdb_from_grid(self, tmpdir):
        """Test basic VDB export from Grid"""
        data = np.arange(1, 28).reshape(3, 3, 3).astype(np.float32)
        g = Grid(data, origin=np.zeros(3), delta=np.ones(3))
        
        outfile = str(tmpdir / "test.vdb")
        g.export(outfile, file_format='VDB')
        
        assert tmpdir.join("test.vdb").exists()
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1
        assert grids[0].name == 'density'

    def test_write_vdb_autodetect_extension(self, tmpdir):
        """Test that .vdb extension is auto-detected"""
        data = np.arange(24).reshape(2, 3, 4).astype(np.float32)
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        
        outfile = str(tmpdir / "auto.vdb")
        g.export(outfile)  
        assert tmpdir.join("auto.vdb").exists()

    def test_write_vdb_default_grid_name(self, tmpdir):
        """Test that default grid name is used when no metadata"""
        data = np.ones((3, 3, 3), dtype=np.float32)
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        g.metadata = {}
        
        outfile = str(tmpdir / "default_name.vdb")
        g.export(outfile)
        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == 'density'

    def test_write_vdb_with_metadata(self, tmpdir):
        """Test that grid name from metadata is used"""
        data = np.ones((3, 3, 3), dtype=np.float32)
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        g.metadata['name'] = 'test_density'
        
        outfile = str(tmpdir / "metadata.vdb")
        g.export(outfile)
        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == 'test_density'

    def test_write_vdb_origin_and_spacing(self, tmpdir):
        """Test that origin and spacing are correctly written"""
        data = np.ones((4, 4, 4), dtype=np.float32)
        origin = np.array([10.0, 20.0, 30.0])
        delta = np.array([0.5, 0.5, 0.5])
        
        g = Grid(data, origin=origin, delta=delta)
        outfile = str(tmpdir / "transform.vdb")
        g.export(outfile)
        
        grids, metadata = vdb.readAll(outfile)
        grid_vdb = grids[0]
        voxel_size = grid_vdb.transform.voxelSize()
        assert_allclose([voxel_size[i] for i in range(3)], delta, rtol=1e-5)

    def test_write_vdb_from_ccp4(self, tmpdir):
        """Test exporting CCP4 file to VDB"""
        g = Grid(datafiles.CCP4)
        outfile = str(tmpdir / "from_ccp4.vdb")
        
        g.export(outfile, file_format='VDB')
        
        assert tmpdir.join("from_ccp4.vdb").exists()
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1

    def test_write_vdb_non3d_raises(self, tmpdir):
        """Test that non-3D grids raise ValueError"""
        data_2d = np.arange(12).reshape(3, 4)
        g = Grid(data_2d, origin=[0, 0], delta=[1, 1])
        
        outfile = str(tmpdir / "invalid.vdb")
        with pytest.raises(ValueError, match="3D grid"):
            g.export(outfile, file_format='VDB')

    def test_write_vdb_with_delta_matrix(self, tmpdir):
        """Test writing with delta as diagonal matrix"""
        data = np.ones((3, 3, 3), dtype=np.float32)
        delta = np.diag([1.0, 2.0, 3.0])
        
        vdb_field = gridData.OpenVDB.field('matrix_delta')
        vdb_field.populate(data, origin=[0, 0, 0], delta=delta)
        
        outfile = str(tmpdir / "matrix_delta.vdb")
        vdb_field.write(outfile)
        assert tmpdir.join("matrix_delta.vdb").exists()

    def test_write_vdb_nonuniform_spacing(self, tmpdir):
        """Test writing with non-uniform spacing"""
        data = np.ones((3, 3, 3), dtype=np.float32)
        delta = np.array([0.5, 1.0, 1.5])
        g = Grid(data, origin=[0, 0, 0], delta=delta)
        
        outfile = str(tmpdir / "nonuniform.vdb")
        g.export(outfile)
        assert tmpdir.join("nonuniform.vdb").exists()

    def test_write_vdb_sparse_data(self, tmpdir):
        """Test writing sparse data"""
        data = np.zeros((10, 10, 10), dtype=np.float32)
        data[2, 3, 4] = 5.0
        data[7, 8, 9] = 10.0
        
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "sparse.vdb")
        g.export(outfile)
        
        assert tmpdir.join("sparse.vdb").exists()
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1

    def test_write_vdb_negative_values(self, tmpdir):
        """Test writing data with negative values"""
        data = np.array([[[1.0, -2.0, 3.0],
                         [-4.0, 5.0, -6.0],
                         [7.0, -8.0, 9.0]]], dtype=np.float32)
        
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "negative.vdb")
        g.export(outfile)
        assert tmpdir.join("negative.vdb").exists()

    def test_write_vdb_zero_threshold(self, tmpdir):
        """Test that very small values near zero are treated as background"""
        data = np.ones((3, 3, 3), dtype=np.float32) * 1e-11 
        data[1, 1, 1] = 1.0 
        
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        outfile = str(tmpdir / "threshold.vdb")
        g.export(outfile)
        assert tmpdir.join("threshold.vdb").exists()

    def test_vdb_field_direct(self, tmpdir):
        """Test using OpenVDB.field directly"""
        data = np.arange(27).reshape(3, 3, 3).astype(np.float32)
        
        vdb_field = gridData.OpenVDB.field('direct_test')
        vdb_field.populate(data, origin=[0, 0, 0], delta=[1, 1, 1])
        
        outfile = str(tmpdir / "direct.vdb")
        vdb_field.write(outfile)
        
        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == 'direct_test'

    def test_vdb_field_no_data_raises(self, tmpdir):
        """Test that writing without data raises ValueError"""
        vdb_field = gridData.OpenVDB.field('empty')
        outfile = str(tmpdir / "empty.vdb")
        with pytest.raises(ValueError, match="No data to write"):
            vdb_field.write(outfile)

    def test_vdb_field_2d_raises(self):
        """Test that 2D data raises ValueError in populate"""
        data_2d = np.arange(12).reshape(3, 4)
        vdb_field = gridData.OpenVDB.field('test')
        
        with pytest.raises(ValueError, match="3D grids"):
            vdb_field.populate(data_2d, origin=[0, 0], delta=[1, 1])


@pytest.mark.skipif(HAS_OPENVDB, reason="Testing import error handling")
def test_vdb_import_error():
    with pytest.raises(ImportError, match="pyopenvdb is required"):
        gridData.OpenVDB.field('test')