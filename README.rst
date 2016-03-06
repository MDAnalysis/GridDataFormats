============================
 README for gridDataFormats
============================

|build| |cov| |docs| |zenodo| |conda|

The **gridDataFormats** package provides classes to unify reading and
writing n-dimensional datasets. One can read grid data from files,
make them available as a `Grid`_ object, and write out the data again.

Availability
------------

The package is licensed under the LGPL, v3 (see files COPYING_ and
`COPYING.LESSER`_) and is available from

* the Python Package Index under the name `GridDataFormats`_
* the GitHub repository https://github.com/MDAnalysis/GridDataFormats

.. _GridDataFormats:
   https://pypi.python.org/pypi/GridDataFormats
.. _COPYING:
   https://raw.githubusercontent.com/MDAnalysis/GridDataFormats/master/COPYING
.. _`COPYING.LESSER`:
   https://raw.githubusercontent.com/MDAnalysis/GridDataFormats/master/COPYING.LESSER
.. _Grid:
   http://www.mdanalysis.org/GridDataFormats/gridData/core.html#gridData.core.Grid

Installation
------------

Install with `pip`_::

  pip install gridDataFormats

.. _pip: https://pip.pypa.io/


Documentation
-------------

For the release docs see the `GridDataFormats docs`_. (The latest docs
are also always available at `griddataformats.readthedocs.org`_.)


.. _GridDataFormats docs:
   http://mdanalysis.org/GridDataFormats
.. _`griddataformats.readthedocs.org`:
   http://griddataformats.readthedocs.org


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

.. |cov| image:: https://coveralls.io/repos/MDAnalysis/GridDataFormats/badge.svg?branch=master&service=github
     :alt: Coverage Status
     :target: https://coveralls.io/github/MDAnalysis/GridDataFormats?branch=master

.. |docs| image:: https://readthedocs.org/projects/griddataformats/badge/?version=latest
    :alt: Documentation
    :target: http://griddataformats.readthedocs.org/en/latest/

.. |zenodo| image:: https://zenodo.org/badge/13219/MDAnalysis/GridDataFormats.svg
    :alt: Zenodo DOI
    :target: https://zenodo.org/badge/latestdoi/13219/MDAnalysis/GridDataFormats

.. |conda| image:: https://anaconda.org/mdanalysis/griddataformats/badges/version.svg
    :alt: Anaconda
    :target: https://anaconda.org/MDAnalysis/griddataformats
