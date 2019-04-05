from __future__ import absolute_import, division

import numpy as np
from numpy.testing import (assert_almost_equal,
                           assert_equal)

from gridData import Grid, CCP4

from . import datafiles

def test_ccp4():
    g = Grid(datafiles.CCP4)
    POINTS = 192
    assert_equal(g.grid.flat, np.arange(1, POINTS+1))
    assert_equal(g.grid.size, POINTS)
    assert_almost_equal(g.delta, [3./4, .5, 2./3])
    assert_equal(g.origin, np.zeros(3))

def test_byteorder():
    with open(datafiles.CCP4, 'rb') as ccp4file:
        flag = CCP4.CCP4._detect_byteorder(ccp4file)
    assert flag in ("@", "=", "<"), "flag {} is not '<'".format(flag)
