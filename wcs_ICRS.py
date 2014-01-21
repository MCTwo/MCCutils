"""
Author: Karen Ng <karenyng@ucdavis.edu>
License: BSD 3-clause
Purpose:
"""
from __future__ import(division, print_function, unicode_literals)
import numpy as np
import astropy as ap
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import ICRS
from astropy.coordinates.coordsystems import SphericalCoordinatesBase


class wcs_ICRS(ap.coordinates.builtin_systems.ICRS,
               SphericalCoordinatesBase):
    """
    As observers the conversion of pixel coordinates and wcs coordinates
    is a routine task. This denotes a class that inherits from all the
    goodie from the astropy.coordinates.builtin_systems.ICRS class and
    automatically do the conversion on initialization of the coordinates
    to avoid woes.
    """
    def __init__(self, *args, **kwargs):
        """
        Superclass :
            astropy.coordinates.builtin_systems.ICRS
            astropy.coordinates.coordsystems.SphericalCoordinatesBase


        Parameters:
        ==========
        See documentation from superclass

        Methods:
        =======
        set_wcs

        """
        ICRS.__init__(self, *args, **kwargs)
        #print self.ra


    def get_WCS_from_fits(self, fitsname, verbose=False):
        """
        Parameters
        ----------
        fitsname = string that contains the name of the fits file
        """
        hdulist = fits.open(fitsname)
        w = wcs.WCS(hdulist[0].header)
        if verbose:
            w.wcs.print_contents()
        return w


    def set_wcs(self, wcs=None, fitsname=None):
        """set
        Parameters:
        ==========

        """
        if fitsname is not None:
            self.wcs = self.get_WCS_from_fits(fitsname)
        elif wcs is not None:
            # store wcs info
            assert type(wcs) == ap.wcs.wcs.WCS, "print not " + \
                " of right  type,\ntype(wcs) needs to be " + \
                "astropy.wcs.wcs.WCS"
            self.wcs = wcs
        else:
            raise ValueError("input of fitsname and wcs are both " +
                             "None - one of them has to be valid")


        # kluegy way of doing transformation
        # want to generalize this to a bunch of ICRS not just one.....
        [self.pixcoord, junk] = self.wcs.wcs_world2pix(np.array(
            [[super(wcs_ICRS, self).ra.value,
              super(wcs_ICRS, self).dec.value],
             [super(wcs_ICRS, self).ra.value,
              super(wcs_ICRS, self).dec.value]]),1)



