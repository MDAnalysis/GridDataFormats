Basic use
=========

In most cases, only one class is important, the
:class:`~gridData.core.Grid`, so we just load this right away::

  from gridData import Grid


Loading data
------------

From a OpenDX file::

  g = Grid("density.dx")

(See also :ref:`opendx-read-write` for more information, especially
when working with visualization programs such as PyMOL, VMD, or
Chimera.)

From a gOpenMol PLT file::

  g = Grid("density.plt")

From the output of :func:`numpy.histogramdd`::

  import numpy
  r = numpy.random.randn(100,3)
  H, edges = np.histogramdd(r, bins = (5, 8, 4))
  g = Grid(H, edges=edges)

For other ways to load data, see the docs for :class:`~gridData.core.Grid`.


Writing out data
----------------

Some formats support writing data (see
:ref:`supported-file-formats` for more details), using the
:meth:`gridData.core.Grid.export` method::

   g.export("density.dx")

The format can also be specified explicitly::

  g.export("density.pkl", file_format="pickle")
   
Some of the exporters (such as for OpenDX, see
:ref:`opendx-read-write`) may take additional, format-specific
keywords, which are documented separately.
	


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

