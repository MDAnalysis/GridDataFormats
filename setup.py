# setuptools installation of Hop
# Copyright (c) 2007-2015 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
from setuptools import setup, find_packages
import versioneer
import codecs

with codecs.open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

setup(name="GridDataFormats",
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="Reading and writing of data on regular grids in Python",
      long_description=long_description,
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="LGPLv3",
      url="https://github.com/MDAnalysis/GridDataFormats",
      download_url="https://github.com/MDAnalysis/GridDataFormats/releases",
      keywords="science array density",
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Microsoft :: Windows ',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      packages=find_packages(exclude=[]),
      package_data={'gridData': ['tests/datafiles/*.dx', 'tests/datafiles/*.ccp4',
                                 'tests/datafiles/*.plt']},
      install_requires=['numpy>=1.0.3', 'six', 'scipy'],
      tests_require=['pytest', 'numpy'],
      zip_safe=True,
      )
