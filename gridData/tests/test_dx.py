from __future__ import absolute_import, division

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

import pytest

import gridData.OpenDX
from gridData import Grid

from . import datafiles

def test_read_dx():
    g = Grid(datafiles.DX)
    POINTS = 8
    ref = np.ones(POINTS)
    ref[4] = 1e-6
    ref[5] = -1e+6
    assert_equal(g.grid.flat, ref)
    assert_equal(g.grid.size, POINTS)
    assert_equal(g.delta, np.ones(3))
    assert_equal(g.origin, np.array([20.1, 3., -10.]))

def test_read_dxgz():
    g = Grid(datafiles.DXGZ, file_format='DXGZ')
    POINTS = 8
    ref = np.ones(POINTS)
    ref[4] = 1e-6
    ref[5] = -1e+6
    assert_equal(g.grid.flat, ref)
    assert_equal(g.grid.size, POINTS)
    assert_equal(g.delta, np.ones(3))
    assert_equal(g.origin, np.array([20.1, 3., -10.]))

@pytest.mark.parametrize("nptype,dxtype", [
    ("float16", "float"),
    ("float32", "float"),
    ("float64", "double"),
    ("int64", "int"),
    ("int32", "int"),
    ("uint32", "unsigned int"),
    ("uint64", "unsigned int"),
    ("int16", "short"),
    ("uint16", "unsigned short"),
    ("int8", "signed byte"),
    ("uint8", "byte"),
])
def test_write_dx(tmpdir, nptype, dxtype, counts=100, ndim=3):
    # conversion from numpy array to DX file

    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    assert_equal(g.grid.sum(), counts)

    with tmpdir.as_cwd():
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

@pytest.mark.parametrize("nptype,dxtype", [
    ("float16", "float"),
    ("float32", "float"),
    ("float64", "double"),
    ("int64", "int"),
    ("int32", "int"),
    ("uint32", "unsigned int"),
    ("uint64", "unsigned int"),
    ("int16", "short"),
    ("uint16", "unsigned short"),
    ("int8", "signed byte"),
    ("uint8", "byte"),
])
def test_write_dxgz(tmpdir, nptype, dxtype, counts=100, ndim=3):
    # conversion from numpy array to DXGZ file

    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    assert_equal(g.grid.sum(), counts)

    with tmpdir.as_cwd():
        outfile = "grid"
        g.export(outfile, file_format='DXGZ')
        outfile += ".dx.gz"
        g2 = Grid(outfile, file_format='DXGZ')

        # check that dxtype was written
        dx = gridData.OpenDX.field(0)
        dx.read(outfile, gz=True)
        data = dx.components['data']
        out_dxtype = data.type

    assert_almost_equal(g.grid, g2.grid,
                        err_msg="written grid does not match original")
    assert_almost_equal(
        g.delta, g2.delta,
        decimal=6,
        err_msg="deltas of written grid do not match original")

    assert_equal(out_dxtype, dxtype)

@pytest.mark.parametrize('nptype', ("complex64", "complex128", "bool_"))
@pytest.mark.filterwarnings("ignore:array dtype.name =")
def test_write_dx_ValueError(tmpdir, nptype, counts=100, ndim=3):
    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    with pytest.raises(ValueError):
        with tmpdir.as_cwd():
            outfile = "grid.dx"
            g.export(outfile)

@pytest.mark.parametrize('nptype', ("complex64", "complex128", "bool_"))
@pytest.mark.filterwarnings("ignore:array dtype.name =")
def test_write_dxgz_ValueError(tmpdir, nptype, counts=100, ndim=3):
    h, edges = np.histogramdd(np.random.random((counts, ndim)), bins=10)
    g = Grid(h, edges)

    # hack the grid to be a different dtype
    g.grid = g.grid.astype(nptype)

    with pytest.raises(ValueError):
        with tmpdir.as_cwd():
            outfile = "grid"
            g.export(outfile, file_format='DXGZ')
