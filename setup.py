# setuptools installation of Hop
# Copyright (c) 2007-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(name="GridDataFormats",
      version="0.2.1",
      description="A python library for reading and writing gridded data",
      long_description="""The gridDataFormats package provides classes to 
unify reading and writing n-dimensional datasets. At the moment this simply
means reading and writing OpenDX files, or reading gOpenMol plt files.
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software/#Hop",
      keywords="science",
      packages=find_packages(exclude=[]),
      package_data = {},
      install_requires=['numpy>=1.0.3',
                        'scipy',          # for remapping/interpolation
                        ], 
      zip_safe=True,
)
