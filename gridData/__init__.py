# $Id$
# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser Public License, version 3 or later.
# See the files COPYING and COPYING.LESSER for details.

"""
:mod:`gridDataFormat` -- Handling grids of data
===============================================

This module contains classes that allow importing and exporting of
simple gridded data, A grid is an N-dimensional array that represents
a discrete mesh over a region of space. The array axes are taken to be
parallel to the cartesian axes of this space. Together with this array
we also store the edges, which are are (essentially) the cartesian
coordinates of the intersections of the grid (mesh) lines on the
axes. In this way the grid is anchored in space.

Reading grid data files
-----------------------

Some Formats_ can be read directly from a file on disk::

 g = Grid(filename)

*filename* could be, for instance, "density.dx".


Constructing a Grid
-------------------

Data from an n-dimensional array can be packaged as a :class:`Grid`
for convenient handling (especially export to other formats).  The
:class:`Grid` class acts as a universal constructor::

 g = Grid(ndarray, edges=edges)                 # from histogramdd
 g = Grid(ndarray, origin=origin, delta=delta)  # from arbitrary data

 g.export(filename, format)   # export to the desire format

See the doc string for :class:`Grid` for details.


Formats
-------

The following formats are available:

   :mod:`OpenDX`
        IBM's Data Explorer, http://www.opendx.org/
   :mod:`gOpenMol`
        http://www.csc.fi/gopenmol/  ## not implemented yet
   pickle
        python pickle file (:mod:`pickle`)

Exceptions
----------

.. autoexception:`gridDataWarning`

"""

__all__ =  ['Grid', 'OpenDX','gOpenMol']

import warnings

class gridDataWarning(Warning):
    """Warns of a problem specific to the gridData module."""
    pass

from core import Grid
