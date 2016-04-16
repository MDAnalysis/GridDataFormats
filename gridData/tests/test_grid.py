import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, dec
from nose.tools import raises

from gridData import Grid

from gridData.tests import module_not_found


class TestGrid:

    def __init__(self):
        self.griddata = np.arange(1, 28).reshape(3, 3, 3)
        self.origin = np.zeros(3)
        self.delta = np.ones(3)
        self.grid = Grid(self.griddata, origin=self.origin, delta=self.delta)

    def test_init(self):
        g = Grid(self.griddata, origin=self.origin, delta=1)
        assert_array_equal(g.delta, self.delta)

    @raises(TypeError)
    def test_init_wrong_origin(self):
        Grid(self.griddata, origin=np.ones(4), delta=self.delta)

    @raises(TypeError)
    def test_init_wrong_delta(self):
        Grid(self.griddata, origin=self.origin, delta=np.ones(4))

    def test_equality(self):
        assert self.grid == self.grid
        assert self.grid != 'foo'
        g = Grid(self.griddata, origin=self.origin + 1, delta=self.delta)
        assert self.grid != g

    def test_addition(self):
        g = self.grid + self.grid
        assert_array_equal(g.grid.flat, (2 * self.griddata).flat)
        g = 2 + self.grid
        assert_array_equal(g.grid.flat, (2 + self.griddata).flat)
        g = g + self.grid
        assert_array_equal(g.grid.flat, (2 + (2 * self.griddata)).flat)

    def test_substraction(self):
        g = self.grid - self.grid
        assert_array_equal(g.grid.flat, np.zeros(27))
        g = 2 - self.grid
        assert_array_equal(g.grid.flat, (2 - self.griddata).flat)

    def test_multiplication(self):
        g = self.grid * self.grid
        assert_array_equal(g.grid.flat, (self.griddata ** 2).flat)
        g = 2 * self.grid
        assert_array_equal(g.grid.flat, (2 * self.griddata).flat)

    def test_division(self):
        # __truediv__ is used in py3 by default and py2 if division
        # is imported from __future__ to make testing easier lets call
        # them explicitely
        g = self.grid.__truediv__(self.grid)
        assert_array_equal(g.grid.flat, np.ones(27))
        g = self.grid.__rtruediv__(2)
        assert_array_equal(g.grid.flat, (2 / self.griddata).flat)

    def test_old_division(self):
        # this is normally ONLY invoked in python 2. To have test
        # coverage in python3 as well call it explicitely
        g = self.grid.__div__(self.grid)
        assert_array_equal(g.grid.flat, np.ones(27))
        g = self.grid.__rdiv__(2)
        assert_array_equal(g.grid.flat, (2 / self.griddata).flat)

    def test_power(self):
        g = self.grid ** 2
        assert_array_equal(g.grid.flat, (self.griddata ** 2).flat)
        g = 2 ** self.grid
        assert_array_equal(g.grid.flat, (2 ** self.griddata).flat)

    def test_compatibility_type(self):
        assert self.grid.check_compatible(self.grid)
        assert self.grid.check_compatible(3)
        g = Grid(self.griddata, origin=self.origin - 1, delta=self.delta)
        assert self.grid.check_compatible(g)

    @raises(TypeError)
    def test_wrong_compatibile_type(self):
        self.grid.check_compatible("foo")

    @raises(NotImplementedError)
    def test_non_orthonormal_boxes(self):
        delta = np.eye(3)
        Grid(self.griddata, origin=self.origin, delta=delta)

    def test_centers(self):
        # this only checks the edges. If you know an alternative
        # algorithm that isn't an exact duplicate of the one in
        # g.centers to test this please implement it.
        g = Grid(self.griddata, origin=np.ones(3), delta=self.delta)
        centers = np.array(list(g.centers()))
        assert_array_equal(centers[0], g.origin)
        assert_array_equal(centers[-1] - g.origin,
                           (np.array(g.grid.shape) - 1) * self.delta)

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def test_resample_factor(self):
        g = self.grid.resample_factor(2)
        assert_array_equal(g.delta, np.ones(3) * .5)
        assert_array_equal(g.grid.shape, np.ones(3) * 6)
        # check that the edges are the same
        assert_array_almost_equal(g.grid[::5, ::5, ::5],
                                  self.grid.grid[::2, ::2, ::2])
