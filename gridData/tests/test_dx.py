import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_allclose

import pytest

import gridData.OpenDX
from gridData import Grid

from . import datafiles

@pytest.mark.parametrize("infile", [datafiles.DX, datafiles.DXGZ])
def test_read_dx(infile):
    g = Grid(infile)
    POINTS = 8
    ref = np.ones(POINTS)
    ref[4] = 1e-6
    ref[5] = -1e+6
    assert_equal(g.grid.flat, ref)
    assert_equal(g.grid.size, POINTS)
    assert_equal(g.delta, np.ones(3))
    assert_equal(g.origin, np.array([20.1, 3., -10.]))

@pytest.mark.parametrize("outfile", ["grid.dx", "grid.dx.gz"])
@pytest.mark.parametrize("nptype,dxtype", [
    ("float16", "float"),
    ("float32", "float"),
    ("float64", "double"),
    ("int64", "int"),
    ("int32", "int"),
    ("uint32", "unsigned int"),
    ("uint64", "unsigned int"),
    ("int16", "short"),
    ("uint16", "unsigned short"),
    ("int8", "signed byte"),
    ("uint8", "byte"),
])
def test_write_dx(tmpdir, nptype, dxtype, outfile, counts=100, ndim=3):
    # conversion from numpy array to DX file

    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    assert_equal(g.grid.sum(), counts)

    with tmpdir.as_cwd():
        g.export(outfile)
        g2 = Grid(outfile)

        # check that dxtype was written
        dx = gridData.OpenDX.field(0)
        dx.read(outfile)
        data = dx.components['data']
        out_dxtype = data.type

    assert_almost_equal(g.grid, g2.grid,
                        err_msg="written grid does not match original")
    assert_almost_equal(
        g.delta, g2.delta,
        decimal=6,
        err_msg="deltas of written grid do not match original")

    assert_equal(out_dxtype, dxtype)

@pytest.mark.parametrize("outfile", ["grid.dx", "grid.dx.gz"])
@pytest.mark.parametrize('nptype', ("complex64", "complex128", "bool_"))
@pytest.mark.filterwarnings("ignore:array dtype.name =")
def test_write_dx_ValueError(tmpdir, nptype, outfile, counts=100, ndim=3):
    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    with pytest.raises(ValueError):
        with tmpdir.as_cwd():
            g.export(outfile)


def test_delta_precision(tmpdir):
    '''Test if the delta has been written to the 7th significant figure.'''
    g = Grid(datafiles.DX)
    g.delta = np.array([90, 90, 150]) / 257
    with tmpdir.as_cwd():
        g.export("grid.dx")
        g2 = Grid("grid.dx")
    assert_almost_equal(
        g.delta, g2.delta,
        decimal=7,
        err_msg="deltas of written grid do not match original")

def test_dx_from_grid():
    data = np.arange(60, dtype=np.float32).reshape((3, 4, 5))
    g = Grid(data, origin=[0, 0, 0], delta=[1, 3, 2.5])
    
    g.metadata['name'] = 'test_density'
    g.metadata['author'] = 'test_user'
    
    dx_field = gridData.OpenDX.field.from_grid(g)
    
    assert isinstance(dx_field, gridData.OpenDX.field)
    
    assert any(('test_density' and 'test_user') in str(c) for c in dx_field.comments) 
    
    assert_allclose(dx_field.components['data'].array, data)
    assert_allclose(dx_field.components['positions'].origin, g.origin)
    
def test_dx_native():
    data = np.arange(60, dtype=np.float32).reshape((3, 4, 5))
    g = Grid(data, origin=[0, 0, 0], delta=[1, 3, 2.5])
    
    dx_field = gridData.OpenDX.field.from_grid(g)
    assert dx_field.native is dx_field
    assert isinstance(dx_field.native, gridData.OpenDX.field)