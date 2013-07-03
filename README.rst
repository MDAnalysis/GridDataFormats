============================
 README for gridDataFormats
============================

The package reads grid data from files, makes them available as a
:class:`Grid` object, and allows one to write out the data again.

A :class:`Grid` consists of a rectangular, regular, N-dimensional
array of data. It contains
(1) The position of the array cell edges.
(2) The array data itself.

This is equivalent to knowing
(1) The origin of the coordinate system (i.e. which data cell
    corresponds to (0,0,...,0)
(2) The spacing of the grid in each dimension.
(3) The data on a grid.

:class:`Grid` objects have some convenient properties:

* The data is represented as a :class:`numpy.array` and thus shares
  all the advantages coming with this sophisticated and powerful
  library.

* They can be manipulated arithmetically, e.g. one can simply add or
  subtract two of them and get another one, or multiply by a
  constant. Note that all operations are defined point-wise (see the
  :mod:`numpy` documentation for details) and that only grids defined
  on the same cell edges can be combined.

* A :class:`Grid` object can also be created from within python code
  e.g. from the output of the :func:`numpy.histogramdd` function.

* The representation of the data is abstracted from the format that
  the files are saved in. This makes it straightforward to add
  additional readers for new formats.

* The data can be written out again in formats that are understood by
  other programs such as VMD or PyMOL.


Examples
========

In most cases, only one class is important, the
:class:`gridData.Grid`, so we just load this right away::

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
  H, edges = numpy.histogramdd(r, bins = (5, 8, 4))
  g = Grid(H, edges=edges)

For other ways to load data, see the docs for :class:`gridData.Grid`.



Subtracting two densities
-------------------------

Assuming one has two densities that were generated on the same grid
positions, stored in files ``A.dx`` and ``B.dx``, one first reads the
data into two :class:`Grid` objects::

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


Classes
=======

.. class:: Grid

   Class to manage a multidimensional grid object.

   A grid is a n-dimensional array and data indicating the dimensions
   of the edges, i.e. the positions of all grid cells in cartesian
   coordinates. 

   :class:`Grid` objects can be multiplied by scalars in order to
   scale the values (but not the dimensions) and like :class:`Grid`
   objects (i.e. corresponding grid points are located at the same
   coordinates) can be used in basic arithmetic operations (addition,
   subtraction, multiplication, division, exponentiation). Basic
   arithmetic is left and right associative, i.e. ``Grid1 + Grid2 ==
   Grid2 + Grid1`` and also ``3 * Grid`` and ``Grid/0.5`` work.

   .. method:: load(filename)

      Load saved (pickled or dx) grid and edges from *filename*.

   .. method:: export(filename[, format])

      Export to file using the given *format*. 

      The format can also be deduced from the suffix of the filename
      though the *format* keyword takes precedence.

      The default format 'dx'.
        
      Implemented formats:

       ========== =================================================
       suffix     description
       ========== =================================================
        dx        OpenDX

        pickle    Python pickle (use :meth:`load` to restore); 
                  :meth:`save` is a short cut for 
                  ``export(format='python')``
       ========== =================================================

   .. method:: save(filename)

      Save a grid object to *filename*.pickle.


   .. method:: resample(edges)

      Resample data to a new grid with edges *edges*.  The order of
      the interpolation is set by
      :attr:`interpolation_spline_order`. Returns a new :class:`Grid`.

   .. attribute:: interpolation_spline_order

      Order of the B-spline interpolation of the data. 3 = cubic; 4 &
      5 are also supported. Only choose values that are acceptable to
      :func:`scipy.ndimage.spline_filter`!

   .. method:: resample_factor(factor)

      Resample to a new regular grid with *factor* * oldN cells along each
      dimension. Returns a new :class:`Grid`.

   .. method:: centers()

      Returns the coordinates of the centers of all grid cells as an
      iterator.

   .. method:: check_compatible(other)

        Check if *other* can be used in an arithmetic operation.

        1) *other* is a scalar
        2) *other* is a grid defined on the same edges
        
        Raises :exc:`TypeError` if not compatible.

   
