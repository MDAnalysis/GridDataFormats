# setuptools installation of Hop
# Copyright (c) 2007-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

long_description = """\

The *gridDataFormats* package provides classes to unify reading and
writing n-dimensional datasets. One can read grid data from files,
make them available as a :class:`Grid` object, and allows one to
write out the data again.

The Grid class
--------------

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
  other programs such as VMD_ or PyMOL_.

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMOL: http://www.pymol.org/


Supported file formats
----------------------

The package can be easily extended. The OpenDX format is widely
understood by many molecular viewers and is sufficient for many
applications that were encountered so far. Hence, at the moment only a
small number of file formats is directly supported.

==========  =========  =====  =====  =========================================
format      extension  read   write  remarks
==========  =========  =====  =====  =========================================
OpenDX_     dx         x      x      subset of OpenDX implemented
gOpenMol    plt        x            
pickle      pickle     x      x      standard Python pickle of the Grid class
==========  =========  =====  =====  =========================================
"""

setup(name="GridDataFormats",
      version="0.2.2",
      description="Reading and writing of data on regular grids in Python",
      long_description=long_description,
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="https://github.com/orbeckst/GridDataFormats",
      keywords="science array density",
      classifiers = ['Development Status :: 4 - Beta', 
                     'Environment :: Console',
                     'Intended Audience :: Science/Research',
                     'License :: OSI Approved :: GNU General Public License (GPL)',
                     'Programming Language :: Python',
                     'Topic :: Scientific/Engineering',
                     'Topic :: Software Development :: Libraries :: Python Modules',
                     ],
      packages=find_packages(exclude=[]),
      package_data = {},
      install_requires=['numpy>=1.0.3',
                        'scipy',          # for remapping/interpolation
                        ], 
      zip_safe=True,
)
