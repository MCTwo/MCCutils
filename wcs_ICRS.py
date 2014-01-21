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
from astropy import units as u
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
        pixcoord = numpy array of numbers with size 2
            denotes the pixcoord
        fitsname = string
            denoting path to fits file containing wcs info
        wcs = astropy.wcs.wcs object
            denote the wcs info


        Methods:
        =======
        set_wcs
        """
        # kluegy way of incorporting optional parameters
        if kwargs.has_key('pixcoord'):
            self.pixcoord = kwargs.get('pixcoord', None)
            del kwargs['pixcoord']

        if kwargs.has_key('fitsname'):
            self.fitsname = kwargs.get('fitsname', None)
            del kwargs['fitsname']

        if kwargs.has_key('wcs'):
            temp = kwargs.get('wcs', None)
            assert type(temp) == ap.wcs.wcs.WCS, "print not " + \
                " of right  type,\ntype(wcs) needs to be " + \
                "astropy.wcs.wcs.WCS"
            self.wcs = temp
            del kwargs['wcs']

        self.verbose = False
        if kwargs.has_key('verbose'):
            self.verbose = kwargs.get('verbose', False)
            del kwargs['verbose']

        if hasattr(self, 'wcs'):
            self.set_wcs(local_wcs=self.wcs, verbose=self.verbose)
        else:
            if hasattr(self, 'fitsname'):
                self.set_wcs(fitsname=self.fitsname,
                                verbose=self.verbose)

        if hasattr(self, 'pixcoord'):
            self.convert_pix2world(verbose=self.verbose)
            kwargs['ra'] = self.wcs_coord[0]
            kwargs['dec'] = self.wcs_coord[1]
            kwargs['unit'] = (u.degree, u.degree)
            ICRS.__init__(self, *args, **kwargs)
        else:
            # kluegy way of incorporting optional parameters
            ICRS.__init__(self, *args, **kwargs)
            if hasattr(self, 'wcs'):
                self.convert_world2pix(verbose=self.verbose)



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

        """
        if local_wcs == None:
            if fitsname is not None:
                self.wcs = self.get_WCS_from_fits(fitsname, verbose)
        if (local_wcs is None and fitsname is None):
            raise ValueError("input of fitsname and wcs are both " +
                            "None - one of them has to be valid")

    def convert_world2pix(self, verbose=False):
        # kluegy way of doing transformation
        # want to generalize this to a bunch of ICRS not just one.....
        [self.pixcoord, junk] = self.wcs.wcs_world2pix(np.array(
            [[super(wcs_ICRS, self).ra.degree,
              super(wcs_ICRS, self).dec.degree],
             [super(wcs_ICRS, self).ra.degree,
              super(wcs_ICRS, self).dec.degree]]),1)
        if verbose:
            print "set self.pixcoord to {0}".format(self.pixcoord)


    def convert_pix2world(self, verbose=False):
        assert self.pixcoord.shape[0] == 2, "wrong length for pixcoord"
        [self.wcs_coord, junk] = self.wcs.wcs_pix2world([self.pixcoord,
                                                         self.pixcoord], 1)
