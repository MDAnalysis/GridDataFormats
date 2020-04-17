from __future__ import absolute_import, division

import pytest

import numpy as np
from numpy.testing import (assert_almost_equal,
                           assert_equal)

from gridData import Grid, CCP4

from . import datafiles

@pytest.fixture(scope="module")
def g():
    return Grid(datafiles.CCP4)

def test_ccp4(g):
    POINTS = 192
    assert_equal(g.grid.flat, np.arange(1, POINTS+1))
    assert_equal(g.grid.size, POINTS)
    assert_almost_equal(g.delta, [3./4, .5, 2./3])
    assert_equal(g.origin, np.zeros(3))


@pytest.fixture(scope="module")
def ccp4data():
    return CCP4.CCP4(datafiles.CCP4_1JZV)

@pytest.mark.parametrize('name,value', [
    ('nc', 96),
    ('nr', 76),
    ('ns', 70),
    ('mode', 2),
    ('ncstart', -4),
    ('nrstart', -23),
    ('nsstart', 102),
    ('nx', 84),
    ('ny', 84),
    ('nz', 160),
    ('xlen', 45.79999923706055),
    ('ylen', 45.79999923706055),
    ('zlen', 89.6500015258789),
    ('alpha', 90.0),
    ('beta', 90.0),
    ('gamma', 90.0),
    ('mapc', 2),
    ('mapr', 1),
    ('maps', 3),
    ('amin', -0.9930942058563232),
    ('amax', 9.050403594970703),
    ('amean', -0.0005801090155728161),
    ('ispg', 92),
    ('nsymbt', 640),
    ('lskflg', 0),
    ('bsaflag', '@'),
    ('skwmat', None),
    ('skwtrn', None),
    ('endianness', 'little'),
    ('arms', 0.4034915268421173),
    ('nlabl', 1),
    ('label', ' Map from fft                                                                   '),
])
def test_ccp4_read_header(ccp4data, name, value):
    if type(value) is float:
        assert_almost_equal(ccp4data.header[name], value, decimal=6)
    else:
        assert_equal(ccp4data.header[name], value)


def test_byteorder():
    with open(datafiles.CCP4, 'rb') as ccp4file:
        flag = CCP4.CCP4._detect_byteorder(ccp4file)
    assert flag in ("@", "=", "<"), "flag {} is not '<'".format(flag)
