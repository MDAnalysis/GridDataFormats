import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_equal

from gridData import Grid


def test_dx():
    g = Grid('gridData/tests/test.dx')
    POINTS = 8
    assert_array_equal(g.grid.flat, np.ones(POINTS))
    assert_equal(g.grid.size, POINTS)
    assert_array_equal(g.delta, np.eye(3))
    assert_array_equal(g.origin, np.zeros(3))
