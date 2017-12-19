import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_equal

from nose.tools import raises

import gridData.OpenDX
from gridData import Grid

from gridData.testing import tempdir


def test_read_dx():
    g = Grid('gridData/tests/test.dx')
    POINTS = 8
    assert_array_equal(g.grid.flat, np.ones(POINTS))
    assert_equal(g.grid.size, POINTS)
    assert_array_equal(g.delta, np.ones(3))
    assert_array_equal(g.origin, np.zeros(3))


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

    assert_array_almost_equal(g.grid,
                              g2.grid,
                              err_msg="written grid does not match original")
    assert_array_almost_equal(
        g.delta,
        g2.delta,
        err_msg="deltas of written grid do not match original")

    assert_equal(out_dxtype, dxtype)


def test_write_dx_float(nptype="float32", dxtype="float"):
    return _test_write_dx(nptype="float32", dxtype="float")

def test_write_dx_double(nptype="float64", dxtype="double"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_int(nptype="int64", dxtype="int"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

def test_write_dx_unsigned_byte(nptype="uint8", dxtype="byte"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)

# more: see OpenDX.array.types

@raises(ValueError)
def test_write_dx_ValueError(nptype="float128", dxtype="unknown"):
    return _test_write_dx(nptype=nptype, dxtype=dxtype)
