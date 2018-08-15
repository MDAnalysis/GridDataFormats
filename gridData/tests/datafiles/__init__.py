from __future__ import absolute_import

from pkg_resources import resource_filename

__all__ = ["DX", "CCP4"]

DX = resource_filename(__name__, 'test.dx')
CCP4 = resource_filename(__name__, 'test.ccp4')
