Installation
============

GridDataFormats can be easily installed via the :ref:`conda<conda-install>` or
:ref:`pip<pip-install>` package managers. 

It is a pure-python package but it has a few other packages (namely
scipy) as dependencies that contain compiled code. For ease of
installation we recommend :ref:`conda<conda-install>` but
:ref:`pip<pip-install>` and installation from source are also fully
supported.



.. _conda-install:

Installing GridDataFormats with ``conda``
-----------------------------------------

The `conda`_ package manager installs, runs, and updates whole
environments with all their dependencies.

Installing *GridDataFormats* from the *conda-forge* channel can be
achieved by adding "conda-forge" to your channels with::

    conda config --add channels conda-forge

Once the *conda-forge* channel has been enabled, *GridDataFormats* can
be installed with::

    conda install griddataformats

Any missing dependencies will be automatically downloaded and
installed in the appropriate versions.

You can later update with ::

    conda update griddataformats


.. _conda: https://docs.conda.io/    
    
.. _pip-install:
    
Installing GridDataFormats with ``pip``
---------------------------------------

Install with `pip`_::

  pip install gridDataFormats

and you can later update with ::

  pip install --upgrade gridDataFormats


`pip` also automatically downloads all missing dependencies and will
attempt to compile them if necessary; this step can fail if you do not
have the correct build environment with the necessary compilers
installed. You should then read the pip_ documentation to learn what
is needed or switch to the :ref:`conda installation<conda-install>`.


.. _pip: https://pip.pypa.io/


