# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
# See the files COPYING and COPYING.LESSER for details.

"""
:mod:`gridData` -- Handling grids of data
=========================================

Overview
--------

This module contains classes that allow importing and exporting of
simple gridded data, A grid is an N-dimensional array that represents
a discrete mesh over a region of space. The array axes are taken to be
parallel to the cartesian axes of this space. Together with this array
we also store the edges, which are are (essentially) the cartesian
coordinates of the intersections of the grid (mesh) lines on the
axes. In this way the grid is anchored in space.

The :class:`~gridData.core.Grid` object can be resampled at arbitrary resolution (by
interpolating the data). Standard algebraic operations are defined for
grids on a point-wise basis (same as for :class:`numpy.ndarray`).


Description
-----------

The package reads grid data from files, makes them available as a
:class:`~gridData.core.Grid` object, and allows one to write out the data again.

A :class:`~gridData.core.Grid` consists of a rectangular, regular, N-dimensional
array of data. It contains

(1) The position of the array cell edges.
(2) The array data itself.

This is equivalent to knowing

(1) The origin of the coordinate system (i.e. which data cell
    corresponds to (0,0,...,0)
(2) The spacing of the grid in each dimension.
(3) The data on a grid.

:class:`~gridData.core.Grid` objects have some convenient properties:

* The data is represented as a :class:`numpy.ndarray` and thus shares
  all the advantages coming with this sophisticated and powerful
  library.

* They can be manipulated arithmetically, e.g. one can simply add or
  subtract two of them and get another one, or multiply by a
  constant. Note that all operations are defined point-wise (see the
  :mod:`numpy` documentation for details) and that only grids defined
  on the same cell edges can be combined.

* A :class:`~gridData.core.Grid` object can also be created from within python code
  e.g. from the output of the :func:`numpy.histogramdd` function.

* The representation of the data is abstracted from the format that
  the files are saved in. This makes it straightforward to add
  additional readers for new formats.

* The data can be written out again in formats that are understood by
  other programs such as VMD or PyMOL.


Reading grid data files
-----------------------

Some Formats_ can be read directly from a file on disk::

 g = Grid(filename)

*filename* could be, for instance, "density.dx".


Constructing a Grid
-------------------

Data from an n-dimensional array can be packaged as a :class:`~gridData.core.Grid`
for convenient handling (especially export to other formats).  The
:class:`~gridData.core.Grid` class acts as a universal constructor::

 g = Grid(ndarray, edges=edges)                 # from histogramdd
 g = Grid(ndarray, origin=origin, delta=delta)  # from arbitrary data

 g.export(filename, format)   # export to the desire format

See the doc string for :class:`~gridData.core.Grid` for details.


Formats
-------

The following formats are available (:ref:`supported-file-formats`):

   :mod:`~gridData.OpenDX`
        IBM's Data Explorer, http://www.opendx.org/
   :mod:`~gridData.gOpenMol`
        http://www.csc.fi/gopenmol/
   pickle
        python pickle file (:mod:`pickle`)



Examples
========

In most cases, only one class is important, the
:class:`~gridData.core.Grid`, so we just load this right away::

  from gridData import Grid


Loading data
------------

From a OpenDX file::

  g = Grid("density.dx")

From a gOpenMol PLT file::

  g = Grid("density.plt")

From the output of :func:`numpy.histogramdd`::

  import numpy
  r = numpy.random.randn(100,3)
  H, edges = np.histogramdd(r, bins = (5, 8, 4))
  g = Grid(H, edges=edges)

For other ways to load data, see the docs for :class:`~gridData.core.Grid`.



Subtracting two densities
-------------------------

Assuming one has two densities that were generated on the same grid
positions, stored in files ``A.dx`` and ``B.dx``, one first reads the
data into two :class:`~gridData.core.Grid` objects::

  A = Grid('A.dx')
  B = Grid('B.dx')

Subtract A from B::

  C = B - A

and write out as a dx file::

  C.export('C.dx')

The resulting file ``C.dx`` can be visualized with any OpenDX-capable
viewer, or later read-in again.


Resampling
----------

Load data::

 A = Grid('A.dx')

Interpolate with a cubic spline to twice the sample density::

 A2 = A.resample_factor(2)

Downsample to half of the bins in each dimension::

 Ahalf = A.resample_factor(0.5)

Resample to the grid of another density, B::

 B = Grid('B.dx')
 A_on_B = A.resample(B.edges)

or even simpler ::

 A_on_B = A.resample(B)

.. Note:: The cubic spline generates region with values that did not
   occur in the original data; in particular if the original data's
   lowest value was 0 then the spline interpolation will probably
   produce some values <0 near regions where the density changed
   abruptly.
"""

from .core import Grid
from . import OpenDX
from . import gOpenMol
from . import CCP4
from . import testing

__all__ = ['Grid', 'OpenDX', 'gOpenMol', 'CCP4', 'testing']
__version__ = '0.3.3'
