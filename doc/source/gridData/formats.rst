.. -*- mode: rst; coding: utf-8 -*-

Formats
=======

A limited number of commonly used formats can be read or written. The
formats are particularly suitable to interface with molecular
visualization tools such as VMD_ or PyMOL_.

Adding new formats is not difficult and user-contributed
format reader/writer can be easily integrated---send a `pull request`_.

.. _supported-file-formats:

Supported file formats
----------------------

The package can be easily extended. The OpenDX_ format is widely
understood by many molecular viewers and is sufficient for many
applications that were encountered so far. Hence, at the moment only a
small number of file formats is directly supported.

.. table:: Available file formats in :mod:`gridData`
   
   ==========  =========  =====  =====  =========================================
   format      extension  read   write  remarks
   ==========  =========  =====  =====  =========================================
   OpenDX_     dx         x      x      subset of OpenDX implemented
   gOpenMol_   plt        x
   pickle      pickle     x      x      standard Python pickle of the Grid class
   ==========  =========  =====  =====  =========================================


.. _pull request: https://github.com/MDAnalysis/GridDataFormats/pulls
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMOL: http://www.pymol.org/
.. _OpenDX: http://www.opendx.org/
.. _gOpenMol: http://www.csc.fi/gopenmol/


Format-specific modules
-----------------------

.. toctree::
   :maxdepth: 1

   formats/OpenDX
   formats/gOpenMol
