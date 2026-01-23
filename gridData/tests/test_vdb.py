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
    def test_write_vdb_from_grid(self, tmpdir):
        data = np.arange(1, 28).reshape(3, 3, 3).astype(np.float32)
        g = Grid(data, origin=np.zeros(3), delta=np.ones(3))
        
        outfile = str(tmpdir / "test.vdb")
        g.export(outfile, file_format='VDB')
        
        assert tmpdir.join("test.vdb").exists()
        
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1
        assert grids[0].name == 'density'

    def test_write_vdb_autodetect_extension(self, tmpdir):
        data = np.arange(24).reshape(2, 3, 4).astype(np.float32)
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        
        outfile = str(tmpdir / "auto.vdb")
        g.export(outfile) 
        
        assert tmpdir.join("auto.vdb").exists()

    def test_write_vdb_with_metadata(self, tmpdir):
        data = np.ones((3, 3, 3), dtype=np.float32)
        g = Grid(data, origin=[0, 0, 0], delta=[1, 1, 1])
        g.metadata['name'] = 'test_density'
        
        outfile = str(tmpdir / "metadata.vdb")
        g.export(outfile)
        
        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == 'test_density'

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
        assert_allclose([voxel_size[i] for i in range(3)], delta, rtol=1e-5)

    def test_write_vdb_from_ccp4(self, tmpdir):
        g = Grid(datafiles.CCP4)
        outfile = str(tmpdir / "from_ccp4.vdb")
        
        g.export(outfile, file_format='VDB')
        
        assert tmpdir.join("from_ccp4.vdb").exists()
        grids, metadata = vdb.readAll(outfile)
        assert len(grids) == 1

    def test_vdb_field_direct(self, tmpdir):
        data = np.arange(27).reshape(3, 3, 3).astype(np.float32)
        
        vdb_field = gridData.OpenVDB.field('direct_test')
        vdb_field.populate(data, origin=[0, 0, 0], delta=[1, 1, 1])
        
        outfile = str(tmpdir / "direct.vdb")
        vdb_field.write(outfile)
        
        grids, metadata = vdb.readAll(outfile)
        assert grids[0].name == 'direct_test'

    def test_vdb_field_no_data_raises(self, tmpdir):
        vdb_field = gridData.OpenVDB.field('empty')
        
        outfile = str(tmpdir / "empty.vdb")
        with pytest.raises(ValueError, match="No data to write"):
            vdb_field.write(outfile)

    def test_vdb_field_2d_raises(self):
        data_2d = np.arange(12).reshape(3, 4)
        vdb_field = gridData.OpenVDB.field('test')
        
        with pytest.raises(ValueError, match="3D grids"):
            vdb_field.populate(data_2d, origin=[0, 0], delta=[1, 1])