# setuptools installation of Hop
# Copyright (c) 2007-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(name="Hop",
      version="0.3.1",
      description="Hop analyses solvent dynamics in molecular dynamics trajectories",
      long_description="""\
Hop performs a 'hopping analysis' of molecules in molecular dynamics
(MD) trajectories. Typically, these molecules are water molecules. The
movement of all waters is tracked as they move between hydration
sites. Hydration sites are interpreted as vertices in a graph, while
movement (or 'hops') between them are taken as edges. In this way
hydration is characterized at a coarse grained level as a network of
hops with rate constants and fluxes derived from the MD simulations.\
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software/#Hop",
      keywords="science 'molecular dynamics' analysis hydration water",
      packages=find_packages(exclude=['tests','extras','doc/examples']),
      package_data = {'vmd': ['*.tcl']},
      install_requires=['numpy>=1.0.3',
                        'scipy',
                        'networkx>1.0',
                        'MDAnalysis>=0.6.3',  # or get 0.6.4-dev from svn
                        ],
      dependency_links = [
        "http://code.google.com/p/mdanalysis/downloads/list",
        ],
      extras_require={
          'plotting': ['matplotlib>=0.91.3', # probably already installed
                       'pygraphviz',         # only needed when plotting, not needed for graph building
                       ]
          'heatmap': ['rpy'],                # optional,used for heatmaps
          },
      zip_safe=True,          # vmdcontrol uses pkg_resources to find vmd tcl script
)
