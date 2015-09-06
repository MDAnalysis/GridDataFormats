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


class TestGrid:

    def __init__(self):
        self.griddata = np.arange(1, 28).reshape(3, 3, 3)
        self.origin = np.zeros(3)
        self.delta = np.ones(3)
        self.grid = Grid(self.griddata, origin=self.origin, delta=self.delta)

    def test_addition(self):
        g = self.grid + self.grid
        assert_array_equal(g.grid.flat, (2 * self.griddata).flat)

    def test_substraction(self):
        g = self.grid - self.grid
        assert_array_equal(g.grid.flat, np.zeros(27))

    def test_multiplication(self):
        g = self.grid * self.grid
        assert_array_equal(g.grid.flat, (self.griddata ** 2).flat)

    def test_division(self):
        g = self.grid / self.grid
        assert_array_equal(g.grid.flat, np.ones(27))

    def test_power(self):
        g = self.grid ** 2
        assert_array_equal(g.grid.flat, (self.griddata ** 2).flat)
