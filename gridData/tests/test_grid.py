from __future__ import absolute_import, division

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, dec

import pytest

from gridData import Grid

class TestGrid(object):
    @staticmethod
    @pytest.fixture(scope="class")
    def data():
        d = dict(
            griddata=np.arange(1, 28).reshape(3, 3, 3),
            origin=np.zeros(3),
            delta=np.ones(3))
        d['grid'] = Grid(d['griddata'], origin=d['origin'],
                         delta=d['delta'])
        return d

    def test_init(self, data):
        g = Grid(data['griddata'], origin=data['origin'],
                 delta=1)
        assert_array_equal(g.delta, data['delta'])

    def test_init_wrong_origin(self, data):
        with pytest.raises(TypeError):
            Grid(data['griddata'], origin=np.ones(4), delta=data['delta'])

    def test_init_wrong_delta(self, data):
        with pytest.raises(TypeError):
            Grid(data['griddata'], origin=data['origin'], delta=np.ones(4))

    def test_equality(self, data):
        assert data['grid'] == data['grid']
        assert data['grid'] != 'foo'
        g = Grid(data['griddata'], origin=data['origin'] + 1, delta=data['delta'])
        assert data['grid'] != g

    def test_addition(self, data):
        g = data['grid'] + data['grid']
        assert_array_equal(g.grid.flat, (2 * data['griddata']).flat)
        g = 2 + data['grid']
        assert_array_equal(g.grid.flat, (2 + data['griddata']).flat)
        g = g + data['grid']
        assert_array_equal(g.grid.flat, (2 + (2 * data['griddata'])).flat)

    def test_substraction(self, data):
        g = data['grid'] - data['grid']
        assert_array_equal(g.grid.flat, np.zeros(27))
        g = 2 - data['grid']
        assert_array_equal(g.grid.flat, (2 - data['griddata']).flat)

    def test_multiplication(self, data):
        g = data['grid'] * data['grid']
        assert_array_equal(g.grid.flat, (data['griddata'] ** 2).flat)
        g = 2 * data['grid']
        assert_array_equal(g.grid.flat, (2 * data['griddata']).flat)

    def test_division(self, data):
        # __truediv__ is used in py3 by default and py2 if division
        # is imported from __future__; to make testing easier lets call
        # them explicitely
        g = data['grid'].__truediv__(data['grid'])
        assert_array_equal(g.grid.flat, np.ones(27))
        g = data['grid'].__rtruediv__(2)
        assert_array_equal(g.grid.flat, (2 / data['griddata']).flat)

    def test_old_division(self, data):
        # this is normally ONLY invoked in python 2. To have test
        # coverage in python3 as well call it explicitely
        g = data['grid'].__div__(data['grid'])
        assert_array_equal(g.grid.flat, np.ones(27))
        g = data['grid'].__rdiv__(2)
        assert_array_equal(g.grid.flat, (2 / data['griddata']).flat)

    def test_power(self, data):
        g = data['grid'] ** 2
        assert_array_equal(g.grid.flat, (data['griddata'] ** 2).flat)
        g = 2 ** data['grid']
        assert_array_equal(g.grid.flat, (2 ** data['griddata']).flat)

    def test_compatibility_type(self, data):
        assert data['grid'].check_compatible(data['grid'])
        assert data['grid'].check_compatible(3)
        g = Grid(data['griddata'], origin=data['origin'] - 1, delta=data['delta'])
        assert data['grid'].check_compatible(g)

    def test_wrong_compatibile_type(self, data):
        with pytest.raises(TypeError):
            data['grid'].check_compatible("foo")

    def test_non_orthonormal_boxes(self, data):
        delta = np.eye(3)
        with pytest.raises(NotImplementedError):
            Grid(data['griddata'], origin=data['origin'], delta=delta)

    def test_centers(self, data):
        # this only checks the edges. If you know an alternative
        # algorithm that isn't an exact duplicate of the one in
        # g.centers to test this please implement it.
        g = Grid(data['griddata'], origin=np.ones(3), delta=data['delta'])
        centers = np.array(list(g.centers()))
        assert_array_equal(centers[0], g.origin)
        assert_array_equal(centers[-1] - g.origin,
                           (np.array(g.grid.shape) - 1) * data['delta'])

    def test_resample_factor(self, data):
        pytest.importorskip('scipy')

        g = data['grid'].resample_factor(2)
        assert_array_equal(g.delta, np.ones(3) * .5)
        assert_array_equal(g.grid.shape, np.ones(3) * 6)
        # check that the edges are the same
        assert_array_almost_equal(g.grid[::5, ::5, ::5],
                                  data['grid'].grid[::2, ::2, ::2])
