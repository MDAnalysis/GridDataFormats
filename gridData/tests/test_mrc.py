from __future__ import absolute_import, division

import pytest

import numpy as np
from numpy.testing import (assert_almost_equal,
                           assert_equal)

from gridData import Grid, mrc

from . import datafiles

@pytest.fixture(scope="module")
def g1():
    return Grid(datafiles.CCP4, file_format="MRC")

@pytest.fixture(scope="module")
def g2():
    data = mrc.MRC()
    data.read(datafiles.CCP4)
    grid, edges = data.histogramdd()
    return Grid(grid=grid, edges=edges)

def test_ccp4_Grid(g1):
    _test_ccp4(g1)

def test_ccp4_mrc(g2):
    _test_ccp4(g2)

def _test_ccp4(g):
    POINTS = 192
    assert_equal(g.grid.flat, np.arange(1, POINTS+1))
    assert_equal(g.grid.size, POINTS)
    assert_almost_equal(g.delta, [3./4, .5, 2./3])
    assert_equal(g.origin, np.zeros(3))



@pytest.fixture(scope="module")
def ccp4data():
    return mrc.MRC(datafiles.CCP4_1JZV)

@pytest.mark.parametrize('name,value', [
    # nx, ny, nz are named nc, nr, ns in the CCP4 module
    ('nx', 96),
    ('ny', 76),
    ('nz', 70),
    ('mode', 2),
    ('nxstart', -4),
    ('nystart', -23),
    ('nzstart', 102),
    ('mx', 84),
    ('my', 84),
    ('mz', 160),
    ('cella', np.rec.array((45.8, 45.8, 89.65),
                        dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])),
    ('cellb', np.rec.array((90., 90., 90.),
                           dtype=[('alpha', '<f4'), ('beta', '<f4'), ('gamma', '<f4')])),
    ('mapc', 2),
    ('mapr', 1),
    ('maps', 3),
    ('dmin', -0.9930942058563232),
    ('dmax', 9.050403594970703),
    ('dmean', -0.0005801090155728161),
    ('ispg', 92),
    ('nsymbt', 640),
    ('nversion', 0),
    ('origin', np.rec.array((0., 0., 0.),
                            dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])),
    ('map', b'MAP '),
    ('machst', np.array([68, 65,  0,  0], dtype=np.uint8)),
    ('rms', 0.40349153),
    ('nlabl', 1),
    ('label', np.array([
        b' Map from fft                                                                   ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                ',
        b'                                                                                '],
                    dtype='|S80'))
])
def test_ccp4_read_header(ccp4data, name, value):
    if type(value) is float:
        assert_almost_equal(ccp4data.header[name], value, decimal=6)
    else:
        assert_equal(ccp4data.header[name], value)

def test_triclinic_ValueError():
    with pytest.raises(ValueError,
                       match="Only orthorhombic unitcells are currently "
                       "supported, not"):
        Grid(datafiles.MRC_EMD3001, file_format="MRC")
