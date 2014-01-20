"""
Author: Karen Ng <karenyng@ucdavis.edu>
License: BSD 3-clause
Purpose:
"""
import astropy as ap
from astropy.coordinates import Angle


class wcs_angle(ap.coordinates.angles.Angle):
    """
    As observers the conversion of pixel coordinates and wcs coordinates is
    a routine task. This denotes a class that inherits from all the goodies
    from the astropy.Angle class and automatically do the conversion on
    initialization of the coordinates to avoid woes.
    """
    def __init__(self, angle, unit = None):
        """
        Superclass : astropy.coordinates.angles.Angle

        Parameters:
        ==========
        angle : astropy.coordinates.Angle object


        Methods:
        =======
        set_wcs


        """
        # initialize the properties of Angle!
        Angle(angle, unit)


    def get_WCS_from_fits(self, fitsname, verbose=True):
        """
        Parameters
        ----------
        fitsname = string that contains the name of the fits file
        """
        hdulist = fits.open(fitsname)
        w = wcs.WCS(hdulist[0].header)
        if verbose:
            print "successfully loaded fits wcs info "
            w.wcs.print_contents()
        return w


    def set_wcs(self, fitsname = None ,wcs = None):

        if fitsname != None:
            self.wcs = self.get_WCS_from_fits(fitsname)
        elif wcs != None:
            # store wcs info
            assert type(wcs) == ap.wcs.wcs.WCS, "print not "+\
                " of right  type,\ntype(wcs) needs to be "+\
                "astropy.wcs.wcs.WCS"
            self.wcs = wcs

        print self.degree
        # make use of the Angle class methods to get the angle
        # in degrees
        # then do conversion by wcs_world2pix
        # then store results as a class variable



