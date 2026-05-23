.. -*- mode: rst; coding: utf-8 -*-

API Interoperability
====================

For a number of formats such as MRC or OpenVDB, GridDataFormats uses external
packages to load data into a native data structure. In addition to using these
data structures to provide **file-level interoperability** by reading files in
one format and writing it in another, GridDataFormats is providing since
release 1.2.0 **API level interoperability** (see issue `#160`_). In this way,
users can directly work with *native* objects representing the data instead of
files. This approach is more efficient in workflows and can make code cleaner.

Typically, GridDataFormats does not work directly with *native objects* (such
as a :class:`openvdb.FloatGrid` for OpenVDB or :class:`mrcfile.mrcfile.MrcFile`
for MRC files) but wraps these classes into *adapter classes* (namely,
:class:`gridData.OpenVDB.OpenVDBField` or :class:`gridData.mrc.MRC`). For some
formats, there is no external library available and the GridDataFormats class
is the "native" object (for instance, :class:`gridData.OpenDX.field` for OpenDX
data).

.. _`#160`: https://github.com/MDAnalysis/GridDataFormats/issues/160
.. _`#161`: https://github.com/MDAnalysis/GridDataFormats/issues/161
.. _`#162`: https://github.com/MDAnalysis/GridDataFormats/issues/162

Converting a Grid object to a native object: ``Grid.convert_to()``
------------------------------------------------------------------

.. versionadded:: 1.2.0

1. The class :class:`~gridData.core.Grid` contains the method
   :meth:`~gridData.core.Grid.convert_to` that creates the *native object*.

2. Each *adapter class* :class:`A` that supports the ``convert_to`` API *must*
   implement

   1. the *classmethod* :class:`A.from_grid` with signature
      ``from_grid(grid: Grid, **kwargs) → A`` that will create the *adapter
      class* from a :class:`Grid` (and use any additional optional arguments
      while ignoring any that it cannot process);
   2. the *attribute* :attr:`A.native` that contains the underlying *native
      object*.

      This attribute can be implemented as a property and should be
      considered read-only, i.e., it is not guaranteed that changing this
      object affects the adapter class :class:`A` although it may do so.

3. The :meth:`A.from_grid` method is listed in the
   :meth:`Grid.converter<gridData.core.Grid.converter>` dictionary.
   
      
For example, given a :class:`~gridData.core.Grid` named ``g``, the following
native objects are produced::

  g = gdf.Grid("density.dx")   # -> gdf.Grid
  mrc = g.convert_to("mrc")    # -> mrcfile.mrcinterpreter.MrcInterpreter
  v = g.convert_to("vdb")      # -> openvdb.GridBase (eg a FloatGrid)  

(See issue `#161`_ for additional background.)
  

Creating a Grid from a native object
------------------------------------

.. note:: **Not implemented yet.** See issue `#162`_ for details.

:class:`~gridData.core.Grid` should be able to consume *native objects* in place of files.


For example, a :class:`mrcfile.mrcfile.MrcFile` instance can be used to
instantiate the :class:`Grid` instance::

  mrc = mrcfile.mrcfile.MrcFile("density.ccp4")
  g = Grid(mrc)  

  
