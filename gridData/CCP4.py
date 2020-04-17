# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
#
# Part of the documentation and format specification:
# Copyright Science and Technologies Facilities Council, 2015.

"""
:mod:`CCP4` --- the CCP4 volumetric data format
===============================================

.. versionadded:: 0.3.0

.. _CCP4: http://www.ccp4.ac.uk/html/maplib.html#description

The module provides a simple implementation of a reader for CCP4_
*ccp4* files. CCP4 files are binary files. The :class:`CCP4` reader tries
to guess the endianess of the file, but this can fail (with a
:exc:`TypeError`); you are on your own in this case.

Only the reader is implemented. If you want to write gridded data use a format
that is more standard, such as OpenDX (see :mod:`OpenDX`).


Background
----------

CCP4 format: http://www.ccp4.ac.uk/html/maplib.html#description

Used to be more carefully documented at
http://lsbr.niams.nih.gov/3demc/3demc_maplib.html but currently this is only
accessible through the Google cache
http://webcache.googleusercontent.com/search?q=cache:KRSvXB0S3dsJ:lsbr.niams.nih.gov/3demc/3demc_maplib.html

Grid data CCP4 file format
--------------------------

Copyright Science and Technologies Facilities Council, 2015.

The overall layout of the file is as follows::

    File header (256 longwords)
    Symmetry information
    Map, stored as a 3-dimensional array

The header is organised as 56 words followed by space for ten 80
character text labels as follows::

   1      NC              # of Columns    (fastest changing in map)
   2      NR              # of Rows
   3      NS              # of Sections   (slowest changing in map)
   4      MODE            Data type
                            0 = envelope stored as signed bytes (from
                                -128 lowest to 127 highest)
                            1 = Image     stored as Integer*2
                            2 = Image     stored as Reals
                            3 = Transform stored as Complex Integer*2
                            4 = Transform stored as Complex Reals
                            5 == 0

                            Note: Mode 2 is the normal mode used in
                                  the CCP4 programs. Other modes than 2 and 0
                                  may NOT WORK

   5      NCSTART         Number of first COLUMN  in map
   6      NRSTART         Number of first ROW     in map
   7      NSSTART         Number of first SECTION in map
   8      NX              Number of intervals along X
   9      NY              Number of intervals along Y
  10      NZ              Number of intervals along Z
  11      X length        Cell Dimensions (Angstroms)
  12      Y length                     "
  13      Z length                     "
  14      Alpha           Cell Angles     (Degrees)
  15      Beta                         "
  16      Gamma                        "
  17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
  18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
  19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
  20      AMIN            Minimum density value
  21      AMAX            Maximum density value
  22      AMEAN           Mean    density value    (Average)
  23      ISPG            Space group number
  24      NSYMBT          Number of bytes used for storing symmetry operators
  25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
  26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                          LSKFLG .ne. 0.
  35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                          Skew transformation is from standard orthogonal
                          coordinate frame (as used for atoms) to orthogonal
                          map frame, as

                                  Xo(map) = S * (Xo(atoms) - t)

  38      future use       (some of these are used by the MSUBSX routines
   .          "              in MAPBRICK, MAPCONT and FRODO)
   .          "   (all set to zero by default)
   .          "
  52          "

  53    MAP             Character string 'MAP ' to identify file type
  54    MACHST          Machine stamp indicating the machine type
                          which wrote file
  55      ARMS            Rms deviation of map from mean density
  56      NLABL           Number of labels being used
  57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)

Symmetry records follow - if any - stored as text as in International
Tables, operators separated by ``*`` and grouped into 'lines' of 80
characters (i.e. symmetry operators do not cross the ends of the
80-character 'lines' and the 'lines' do not terminate in a ``*``).

Map data array follows.

Note on the machine stamp: The machine stamp (word 54) is a 32-bit
quantity containing a set of four 'nibbles' (half-bytes) - only half
the space is used. Each nibble is a number specifying the
representation of (in C terms) double (d), float (f), int (i) and
unsigned char (c) types. Thus each stamp is of the form 0xdfic0000.
For little endian hardware the stamp is 0x44, 0x41, 0x00, 0x00 while
the big endian stamp is 0x11, 0x11, 0x00, 0x00.

Classes
-------

"""
from __future__ import absolute_import, division

import warnings
import struct
import numpy as np
from six.moves import range

from .gOpenMol import Record


# TODO: Consider abstracting a binary data class to handle CCP4,
# gOpenMol, and other binary formats.
class CCP4(object):
    """A class to represent a CCP4_ file.

    Only reading is implemented; either supply a filename to the constructor
      >>> G = CCP4(filename)
    or load the file with the read method
      >>> G = CCP4()
      >>> G.read(filename)

    The data is held in :attr:`CCP4.array` and all header information is in
    the dict :attr:`CCP4.header`.

    :attr:`CCP4.shape`
         D-tuplet describing size in each dimension
    :attr:`CCP4.origin`
         coordinates of the centre of the grid cell with index 0,0,...,0
    :attr:`CCP4.delta`
         DxD array describing the deltas

    .. Note:: The following features of the CCP4 format are *not* implemented:
       * triclinic boxes
       * symmetry records
       * index ordering besides standard column-major and row-major
       * non-standard fields, such any in filed in future use block
    """

    _axis_map = {1: 'x', 2: 'y', 3: 'z'}

    _data_bintype = 'f'

    _header_struct = (
        Record('nc', 'I'),  # of columns (fastest varying index.)
        Record('nr', 'I'),  # of rows
        Record('ns', 'I'),  # of sections (slowest varying index.)
        Record('mode', 'I', {
            0: 'envelope',
            1: 'Image of Integer*2',
            2: 'Image of Reals',  # Default expected value.
            3: 'Transform of Complex Integer*2',
            4: 'Transform of Complex Reals',
            5: '0',
        }), Record('ncstart', 'i'), Record('nrstart', 'i'),
        Record('nsstart', 'I'), Record('nx', 'I'),  # Number of gridpoints.
        Record('ny', 'I'), Record('nz', 'I'), Record('xlen', 'f'),  # Angstroms.
        Record('ylen', 'f'), Record('zlen', 'f'), Record('alpha', 'f'),  # Degrees.
        Record('beta', 'f'), Record('gamma', 'f'),
        Record('mapc', 'I', _axis_map), Record('mapr', 'I', _axis_map),
        Record('maps', 'I', _axis_map), Record('amin', 'f'),
        Record('amax', 'f'), Record('amean', 'f'), Record('ispg', 'I'),
        Record('nsymbt', 'I'), Record('lskflg', 'I'),  # Remaining few fields are manually parsed.
    )

    def __init__(self, filename=None):
        self.filename = str(filename)
        # Assemble format.
        self._headerfmt = "".join([r.bintype for r in self._header_struct])

        if filename is not None:
            self.read(filename)

    def read(self, filename):
        """Populate the instance from the ccp4 file *filename*."""
        if filename is not None:
            self.filename = str(filename)
        with open(self.filename, 'rb') as ccp4:
            h = self.header = self._read_header(ccp4)
            nentries = h['nc'] * h['nr'] * h['ns']
            # Quick and dirty... slurp it all in one go.
            datafmt = h['bsaflag'] + str(nentries) + self._data_bintype
            a = np.array(struct.unpack(datafmt, ccp4.read(struct.calcsize(datafmt))))
        self.header['filename'] = self.filename
        # TODO: Account for the possibility that y-axis is fastest or
        # slowest index, which unfortunately is possible in CCP4.
        order = 'C' if h['mapc'] == 'z' else 'F'
        self.array = a.reshape(h['nc'], h['nr'], h['ns'], order=order)
        self.delta = self._delta()
        self.origin = np.zeros(3)
        self.rank = 3

    @property
    def shape(self):
        return self.array.shape

    @property
    def edges(self):
        """Edges of the grid cells, origin at centre of 0,0,..,0 grid cell.

        Only works for regular, orthonormal grids.
        """
        # TODO: Add triclinic cell support.
        return [self.delta[d, d] * np.arange(self.shape[d] + 1) +
                self.origin[d] - 0.5 * self.delta[d, d]
                for d in range(self.rank)]

    def _delta(self):
        h = self.header
        lengths = np.array([h['xlen'], h['ylen'], h['zlen']])
        delta = lengths / self.shape
        return np.diag(delta)

    @staticmethod
    def _detect_byteorder(ccp4file):
        """Detect the byteorder of stream `ccp4file` and return format character.

        Try all endinaness and alignment options until we find
        something that looks sensible ("MAPS " in the first 4 bytes).

        (The ``machst`` field could be used to obtain endianness, but
        it does not specify alignment.)

        .. SeeAlso::

          :mod:`struct`
        """
        bsaflag = None
        ccp4file.seek(52 * 4)
        mapbin = ccp4file.read(4)
        for flag in '@=<>':
            mapstr = struct.unpack(flag + '4s', mapbin)[0].decode('utf-8')
            if mapstr.upper() == 'MAP ':
                bsaflag = flag
                break  # Only possible value according to spec.
        else:
            raise TypeError(
                "Cannot decode header --- corrupted or wrong format?")
        ccp4file.seek(0)
        return bsaflag

    def _read_header(self, ccp4file):
        """Read header bytes"""

        bsaflag = self._detect_byteorder(ccp4file)

        # Parse the top of the header (4-byte words, 1 to 25).
        nheader = struct.calcsize(self._headerfmt)
        names = [r.key for r in self._header_struct]
        bintopheader = ccp4file.read(25 * 4)

        def decode_header(header, bsaflag='@'):
            h = dict(zip(names, struct.unpack(bsaflag + self._headerfmt,
                                              header)))
            h['bsaflag'] = bsaflag
            return h

        header = decode_header(bintopheader, bsaflag)
        for rec in self._header_struct:
            if not rec.is_legal_dict(header):
                warnings.warn(
                    "Key %s: Illegal value %r" % (rec.key, header[rec.key]))

        # Parse the latter half of the header (4-byte words, 26 to 256).
        if (header['lskflg']):
            skewmatrix = np.fromfile(ccp4file, dtype=np.float32, count=9)
            header['skwmat'] = skewmatrix.reshape((3, 3))
            header['skwtrn'] = np.fromfile(ccp4file, dtype=np.float32, count=3)
        else:
            header['skwmat'] = header['skwtrn'] = None
            ccp4file.seek(12 * 4, 1)
        ccp4file.seek(15 * 4, 1)  # Skip future use section.
        ccp4file.seek(4, 1)  # Skip map text, already used above to verify format.
        # TODO: Compare file specified endianness to one obtained above.
        endiancode = struct.unpack(bsaflag + '4b', ccp4file.read(4))
        header['endianness'] = 'little' if endiancode == (0x44, 0x41, 0, 0
                                                          ) else 'big'
        header['arms'] = struct.unpack(bsaflag + 'f', ccp4file.read(4))[0]
        header['nlabl'] = struct.unpack(bsaflag + 'I', ccp4file.read(4))[0]
        if header['nlabl']:
            binlabel = ccp4file.read(80 * header['nlabl'])
            flag = bsaflag + str(80 * header['nlabl']) + 's'
            label = struct.unpack(flag, binlabel)[0]
            header['label'] = label.decode('utf-8').rstrip('\x00')
        else:
            header['label'] = None
        ccp4file.seek(256 * 4)
        # TODO: Parse symmetry records, if any.
        return header

    def histogramdd(self):
        """Return array data as (edges,grid), i.e. a numpy nD histogram."""
        return (self.array, self.edges)
