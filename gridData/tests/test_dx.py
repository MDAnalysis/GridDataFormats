import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_equal

from gridData import Grid
from gridData.testing import tempdir


def test_read_dx():
    g = Grid('gridData/tests/test.dx')
    POINTS = 8
    assert_array_equal(g.grid.flat, np.ones(POINTS))
    assert_equal(g.grid.size, POINTS)
    assert_array_equal(g.delta, np.ones(3))
    assert_array_equal(g.origin, np.zeros(3))


def test_write_dx(counts=100, ndim=3):
    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)
    assert_equal(g.grid.sum(), counts)

    with tempdir.in_tempdir():
        outfile = "grid.dx"
        g.export(outfile)
        g2 = Grid(outfile)

    assert_array_almost_equal(g.grid,
                              g2.grid,
                              err_msg="written grid does not match original")
    assert_array_almost_equal(
        g.delta,
        g2.delta,
        err_msg="deltas of written grid do not match original")
