# gridData --- python modules to read and write gridded data
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Copyright CSC, 2005

"""NOT IMPLEMENTED YET

gOpenMol http://www.csc.fi/english/pages/g0penMol plt format.

Used to be documented at http://www.csc.fi/gopenmol/developers/plt_format.phtml but currently this is only accessible through the internet archive at
http://web.archive.org/web/20061011125817/http://www.csc.fi/gopenmol/developers/plt_format.phtml



Grid data plt file format

Plot file (plt) format The plot files are regular 3D grid files for
plotting of molecular orbitals, electron densities or other molecular
properties. The plot files are produced by several programs. It is
also possible to format/unformat plot files using the pltfile program
in the utility directory. It is also possible to produce plot files
with external (own) programs. Produce first a formatted text file and
use then the pltfile program to unformat the file for gOpenMol. The
format for the plot files are very simple and a description of the
format can be found elsewhere in this manual. gOpenMol can read binary
plot files from different hardware platforms independent of the system
type (little or big endian machines).

Format of the binary *.plt file

The *.plt grid data files are the main enginees for ploting isocontour
surfaces. The *.plt files can either be generated directly by various
programs inside the gOpenMol package or be imported from other
programs.

Programs inside gOpenMol generating *.plt grid data files:

    * The probsurf program generates the grid data for plotting
      Connolly type of surfaces.
    
    * The vss program generates the grid data for plotting
      electrostatic potentials.

To assist in moving the binary *.plt files between different hardware
platforms the pltfile program is included. Using the program it is
possible to format a binary grid file and move it to an other platform
and unformat the file again. The program can of course also be used
interface the grid data from other programs with gOpenMol.

The *.plt file binary and formatted file formats are very simple but
please observe that unformatted files written with a FORTRAN program
are not pure binary files because there are file records between the
values while pure binary files do not have any records between the
values. gOpenMol should be able to figure out if the file is pure
binary or FORTRAN unformatted but it is not very well tested.

Binary *.plt (grid) file format

Record number and meaning:

#1: Integer, rank value must always be = 3
#2: Integer, possible values are 1 ... 50. This value is not used but
it can be used to define the type of surface!
Values used (you can use your own value between 1... 50):

1:   VSS surface
2:   Orbital/density surface
3:   Probe surface
200: Gaussian 94/98
201: Jaguar
202: Gamess
203: AutoDock
204: Delphi/Insight
205: Grid

Value 100 is reserved for grid data coming from OpenMol!

#3: Integer, number of points in z direction
#4: Integer, number of points in y direction
#5: Integer, number of points in x direction
#6: Float, zmin value
#7: Float, zmax value
#8: Float, ymin value
#9: Float, ymax value
#10: Float, xmin value
#11: Float, xmax value
#12 ... Float, grid data values running (x is inner loop, then y and last z):

1.	Loop in the z direction
2.	Loop in the y direction
3.	Loop in the x direction

The formatted (the first few lines) file can look like:

3 2
65 65 65
-3.300000e+001 3.200000e+001 -3.300000e+001 3.200000e+001 -3.300000e+001 3.200000e+001
-1.625609e+001 -1.644741e+001 -1.663923e+001 -1.683115e+001 -1.702274e+001 -1.721340e+001 
-1.740280e+001 -1.759018e+001 -1.777478e+001 -1.795639e+001 -1.813387e+001 -1.830635e+001 
...

Formatted *.plt (grid) file format

Line numbers and variables on the line:
#1: Integer, Integer. Rank and type of surface (rank is always = 3)
#2: Integer, Integer, Integer. Zdim, Ydim, Xdim (number of points in the z,y,x directions)
#3: Float, Float, Float, Float, Float, Float. Zmin, Zmax, Ymin, Ymax, Xmin,Xmax (min and max values)
#4: ... Float. Grid data values running (x is inner loop, then y and last z) with one or several values per line:

1. Loop in the z direction
2. Loop in the y direction
3. Loop in the x direction

A file in this format can then be converted into binary with the
pltfile program.

For an example of what a formatted *.plt file can look like please
look above. If you are a programmer please look at the included
utility programs for the code for doing the reading and writing of
*.plt files

Copyright CSC, 2005
Last modified: September 23, 2003 09:18:50
"""

class GOpenMol(object):
    raise NotImplementedError
