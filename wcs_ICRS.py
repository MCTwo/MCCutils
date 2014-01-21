"""
Author: Karen Ng <karenyng@ucdavis.edu>
License: BSD 3-clause
Purpose:
"""
from __future__ import(division, unicode_literals)
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

        Optional parameters:
        ===================
        fitsname = string
            denoting path to fits file containing wcs info
        wcs = astropy.wcs.wcs object
            denote the wcs info

        Methods:
        =======
        set_wcs
        """
        # kluegy way of incorporting optional parameters
        if kwargs.has_key('fitsname'):
            self.fitsname = kwargs.get('fitsname', None)
            del kwargs['fitsname']

        if kwargs.has_key('wcs'):
            self.wcs = kwargs.get('wcs', None)
            del kwargs['wcs']

        self.verbose = False
        if kwargs.has_key('verbose'):
            self.verbose = kwargs.get('verbose', False)
            del kwargs['verbose']


        ICRS.__init__(self, *args, **kwargs)

        # kluegy way of incorporting optional parameters
        if hasattr(self, 'wcs'):
            self.set_wcs(local_wcs=self.wcs, verbose=self.verbose)
        if hasattr(self, 'fitsname'):
            self.set_wcs(fitsname=self.fitsname, verbose=self.verbose)


    def get_WCS_from_fits(self, fitsname, verbose=False):
        """
        Parameters
        ----------
        fitsname = string that contains the name of the fits file
        """
        hdulist = fits.open(fitsname)
        w = wcs.WCS(hdulist[0].header)
        if verbose:
            print "successfully loaded wcs from fits"
            w.wcs.print_contents()
        return w


    def set_wcs(self, fitsname=None, local_wcs=None, verbose=False):
        """set
        Parameters:
        ==========

        """
        if fitsname is not None:
            self.wcs = self.get_WCS_from_fits(fitsname, verbose)
        else:
            if local_wcs is not None:
                # store wcs info
                assert type(local_wcs) == ap.wcs.wcs.WCS, "print not " + \
                    " of right  type,\ntype(wcs) needs to be " + \
                    "astropy.wcs.wcs.WCS"
                self.wcs = local_wcs
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
        if verbose:
            print "set self.pixcoord to {0}".format(self.pixcoord)
