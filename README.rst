============================
 README for GridDataFormats
============================

|build| |cov| |docs| |zenodo| |conda|

The **GridDataFormats** package provides classes to unify reading and
writing n-dimensional datasets. One can read grid data from files,
make them available as a `Grid`_ object, and write out the data again.

Availability
------------

The package is licensed under the LGPL, v3 (see files COPYING_ and
`COPYING.LESSER`_) and is available 

* from the Python Package Index under the name `GridDataFormats`_
* as a conda package from the *conda-forge* channel, `conda-forge/griddataformats`_ 
* in source from the GitHub repository https://github.com/MDAnalysis/GridDataFormats

.. _GridDataFormats:
   https://pypi.python.org/pypi/GridDataFormats
.. _`conda-forge/griddataformats`:
   https://anaconda.org/conda-forge/griddataformats
.. _COPYING:
   https://raw.githubusercontent.com/MDAnalysis/GridDataFormats/master/COPYING
.. _`COPYING.LESSER`:
   https://raw.githubusercontent.com/MDAnalysis/GridDataFormats/master/COPYING.LESSER
.. _Grid:
   https://www.mdanalysis.org/GridDataFormats/gridData/core.html#gridData.core.Grid

Installation
------------

Installing GridDataFormats with ``pip``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install with `pip`_::

  pip install gridDataFormats

.. _pip: https://pip.pypa.io/


Installing GridDataFormats with ``conda``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing *GridDataFormats* from the *conda-forge* channel can be
achieved by adding "conda-forge" to your channels with::

    conda config --add channels conda-forge

Once the *conda-forge* channel has been enabled, *GridDataFormats* can
be installed with::

    conda install griddataformats



Documentation
-------------

For the latest docs see the `GridDataFormats docs`_. (Multiple
versions of the docs are also available at
`griddataformats.readthedocs.org`_.)


.. _GridDataFormats docs:
   https://www.mdanalysis.org/GridDataFormats
.. _`griddataformats.readthedocs.org`:
   https://griddataformats.readthedocs.org


Contributing
------------

Please use the `issue tracker`_ for bugs and questions.

**GridDataFormats** is open source and contributions are
welcome. Please fork the repository and submit a `pull request`_.

.. _issue tracker:
   https://github.com/MDAnalysis/GridDataFormats/issues
.. _pull request:
   https://github.com/MDAnalysis/GridDataFormats/pulls

.. |build| image:: https://travis-ci.org/MDAnalysis/GridDataFormats.svg?branch=master
    :alt: Build Status
    :target: https://travis-ci.org/MDAnalysis/GridDataFormats

.. |cov| image:: https://codecov.io/gh/MDAnalysis/GridDataFormats/branch/master/graph/badge.svg
     :alt: Coverage Status
     :target: https://codecov.io/gh/MDAnalysis/GridDataFormats

.. |docs| image:: https://readthedocs.org/projects/griddataformats/badge/?version=latest
    :alt: Documentation
    :target: http://griddataformats.readthedocs.org/en/latest/

.. |zenodo| image:: https://zenodo.org/badge/13219/MDAnalysis/GridDataFormats.svg
    :alt: Zenodo DOI
    :target: https://zenodo.org/badge/latestdoi/13219/MDAnalysis/GridDataFormats

.. |conda| image:: https://anaconda.org/conda-forge/griddataformats/badges/version.svg
    :alt: Anaconda
    :target: https://anaconda.org/conda-forge/griddataformats
