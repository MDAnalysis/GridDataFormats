import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

from nose.tools import raises

import gridData.OpenDX
from gridData import Grid

from gridData.testing import tempdir


def test_read_dx():
    g = Grid('gridData/tests/test.dx')
    POINTS = 8
    ref = np.ones(POINTS)
    ref[4] = 1e-6
    ref[5] = -1e+6
    assert_equal(g.grid.flat, ref)
    assert_equal(g.grid.size, POINTS)
    assert_equal(g.delta, np.ones(3))
    assert_equal(g.origin, np.array([20.1, 3., -10.]))


def _test_write_dx(counts=100, ndim=3, nptype="float32", dxtype="float"):
    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    assert_equal(g.grid.sum(), counts)

    with tempdir.in_tempdir():
        outfile = "grid.dx"
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


# conversion from numpy array to DX file

def test_write_dx_float_float16(nptype="float16", dxtype="float"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_float_float32(nptype="float32", dxtype="float"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_double_float64(nptype="float64", dxtype="double"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_int_int64(nptype="int64", dxtype="int"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_int_int32(nptype="int32", dxtype="int"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_unsigned_int_uint32(nptype="uint32", dxtype="unsigned int"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_unsigned_int_uint64(nptype="uint64", dxtype="unsigned int"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_short_int16(nptype="int16", dxtype="short"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_unsigned_short_uint16(nptype="uint16", dxtype="unsigned short"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_signed_byte(nptype="int8", dxtype="signed byte"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_byte(nptype="uint8", dxtype="byte"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)


@raises(ValueError)
def test_write_dx_ValueError(nptype="longdouble", dxtype="unknown"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)
