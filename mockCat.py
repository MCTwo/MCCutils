'''
Author: Karen Y. Ng <karenyng@ucdavis.edu>
Date: 10/05/2013
License: 3 clause BSD
Note:
To be used in conjunction with the functions in nfwMCMC.py by Will Dawson

Dependencies:
py27-astropy @0.3_0
py27-numpy @1.8.0_2
should also note the commit version of homebrewed code just to be safe
and the branch info
'''
from __future__ import division
import numpy
import numpy.random
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy import wcs
import astropy.io.fits
from astropy.coordinates import angle_utilities as ang_util
import pandas as pd

# import homebrewed modules ---------------
from nfwMCMC import shear
#import metro
from profiles import nfwparam
import cosmo
import tools
#import cosmoplot as cplot
#import ellip


def makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd,
                halos, filename, r_bounds=(5e-3, 100), noise_sigma=0.0,
                pix_coord=('x', 'y'), h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0,
                cd=((-1, 0), (0, 1)), verbose=False):

    # read in the catalogs
    cat = tools.readcatalog(catalog)
    key = tools.readheader(catalog)
    # number of galaxies
    N_gal = numpy.shape(cat)[0]

    #-----------Perform some basic user input checks----------------
    # Determine the number of halo's to model
    N_halo = numpy.size(halos) // 3
    if numpy.size(halos) % 3. != 0:
        print 'fitM200: an invalid "halos" input was specified, exiting'
        sys.exit()

    if verbose:
        print 'fitM200: calculating properties of the galaxies' + \
            ' with respect to each halo'

    #-----------Initialize some of the vectors ---------------------
    del_a = numpy.zeros((N_gal, N_halo))
    del_d = numpy.zeros((N_gal, N_halo))
    invSigmacr = numpy.zeros((N_gal, N_halo))
    #var_invSigmacr = numpy.zeros((N_gal,N_halo))
    mask_r_inner = numpy.zeros((N_gal, N_halo)) != 0
    mask_r_outer = numpy.zeros((N_gal, N_halo)) != 0
    mask_z = numpy.zeros((N_gal, N_halo)) != 0

    for h in numpy.arange(N_halo):
        if verbose:
            print 'fitM200: evaluating galaxies with respect' + \
                ' to halo {0}'.format(h)
        for g in numpy.arange(N_gal):
            # Calculate the angular separations
            del_a[g, h], del_d[g, h] = \
                tools.angcomp(cat[g, key[coord[0]]], cat[g, key[coord[1]]],
                              halos[h * 3], halos[h * 3 + 1])

            # *# Note that I could speed up the code by masking galaxies
            # outside
            # *# the apodizing radii and not calculating their invSigmacr
            # Calulate inverse sigma critical
            if halos[h * 3 + 2] < cat[g, key[coord[2]]]:
    #            print 'right before invsigcr calc'
    #            print 'z_halo = {0}'.format(halos[h*3+2])
    #            print 'z_gal = {0}'.format(cat[g,key[coord[2]]])
                invSigmacr[g, h] =\
                    1 / cosmo.lensgeo(halos[h * 3 + 2],
                                      cat[g, key[coord[2]]],
                                      h_scale, Om, Ol)['sigcr']
    #            print 'right after invsigcr calc'
    #            print 'break point'
            # Calculate the standard error on sigma critcal
            #
            # Need to add this section
            #
        # convert deg anglular separation to Mpc separation
        del_a[:, h] *=\
            60 * cosmo.ProjectedLength(halos[h * 3 + 2], h_scale, Om, Ol)
        del_d[:, h] *=\
            60 * cosmo.ProjectedLength(halos[h * 3 + 2], h_scale, Om, Ol)
    del_r = numpy.sqrt(del_a ** 2 + del_d ** 2)
    for h in numpy.arange(N_halo):
        # Mask galaxies that are outside the apodizing radii or
        # foreground of the halo
        mask_r_inner[:, h] = del_r[:, h] > r_bounds[0]
        # mask out galaxies that are not within any outer radii of a halo
        mask_r_outer[:, h] = del_r[:, h] < r_bounds[1]
    #    mask_r[:,h] = numpy.logical_and(del_r[:,h]>r_bounds[0],
    #                                    del_r[:,h]<r_bounds[1])
        # this masks the galaxies that are in front of the
        mask_z[:, h] = cat[:, key[coord[2]]] > halos[h * 3 + 2]
    # mask out all the galaxies that are within any of the inner radii
    # of either halo
    mask_r_inner = numpy.sum(mask_r_inner, axis=1) == N_halo

    #mask = mask_r*mask_z
    mask = numpy.reshape(mask_r_inner, (numpy.size(mask_r_inner), 1)) *\
        mask_r_outer * mask_z

    #mask = mask_r*mask_z
    mask = numpy.reshape(mask_r_inner, (numpy.size(mask_r_inner), 1)) *\
        mask_r_outer * mask_z

    if verbose:
        for h in numpy.arange(N_halo):
            print 'fitM200: There are {0}'.format(numpy.sum(mask[:, h])) +\
                ' background galaxies for Halo {1}'.format(h)

    # there are some initialization and calculation missing
    # here !!!!!!

    # Find the unit relation between del_x, del_y and del_a, del_d

    cd_inv = numpy.linalg.inv(numpy.array(cd))
    cd_mag = numpy.sqrt(numpy.sum(cd_inv[0, :] ** 2))
    del_x = cd_inv[0, 0] / cd_mag * del_a + cd_inv[0, 1] / cd_mag * del_d
    del_y = cd_inv[1, 0] / cd_mag * del_a + cd_inv[1, 1] / cd_mag * del_d
    #------Calculate the  angle of the galaxy CCW from +x axis of the halo
    #-------use trig identities to ---------------------
    #-------sinphi = del_y/del_r -----------------------
    #-------cos2phi = 1-2sinphi**2 ---------------------
    #-------sin2phi = 2sinphi*cosphi -------------------
    cos2phi = 1 - 2 * (del_y / del_r) ** 2
    sin2phi = 2 * del_x * del_y / del_r ** 2

    N_halo = numpy.shape(del_r)[1]
    N_gal = numpy.shape(del_r)[0]
    gamma = numpy.zeros((N_gal, N_halo))
    z_halo = numpy.zeros(N_halo)
    del_c = numpy.zeros(N_halo)
    r_s = numpy.zeros(N_halo)
    for h in numpy.arange(N_halo):
        z_halo[h] = halos[h * 3 + 2]
        # print z_halo[h]

    # Calculate the NFW parameters corresponding to M_200
    # print 'nfwMCMC.nfwparam - M_200 is ',M_200
    for h in numpy.arange(N_halo):
        del_c[h], r_s[h] = nfwparam(M_200_pred[h], z_halo[h], h_scale=0.7,
                                    Om=0.3, Ol=0.7, Or=0.0)
        # loop through each halo calculating the expected absolute shear
        # due to that individual halo
        # *# Could place a mapping based multiprocessing function here
        # calculate the expected ellipticity components of each galaxy
        gamma[mask[:, h], h] = \
            shear(del_c[h], r_s[h], del_r[mask[:, h], h],
                  z_halo[h], invSigmacr[mask[:, h], h],
                  h_scale, Om, Ol, Or)
    mask_e = numpy.sum(mask, axis=1) != 0
    #--calculate the expected ellipticities of each galaxy due to all halos
    e1_exp = numpy.sum(-gamma[mask_e, :] * cos2phi[mask_e, :], axis=1)
    e2_exp = numpy.sum(-gamma[mask_e, :] * sin2phi[mask_e, :], axis=1)

    e1_exp = addNoise(e1_exp, noise_sigma)
    e2_exp = addNoise(e2_exp, noise_sigma)

    # need to check the dimensionality of e1_exp and e2_exp after applying
    # the masks
    #---write out the ellipiticies and the corresponding RA / DEC / X / Y
    #---coordinate to a ttype header list
    headernames = (coord[0], coord[1], coord[2], ellip_theo[0],
                   ellip_theo[1], ellip_theo[2])
    # this is the way to call a column in the catalog e.g.
    # cat[:,key['name']

    #----------figure out all the key names --------
    RA_col = key[coord[0]]
    DEC_col = key[coord[1]]
    z_col = key[coord[2]]
    #    e1_col = ellip_theo[0]
    #    e2_col = ellip_theo[1]
    #    de_col = ellip_theo[2]

    #-------initialize an array that denotes std dev of the ellipticities
    de_exp = np.zeros(cat.shape[0])

    # stack all the catalog content together in one big 2D numpy array
    # remember to apply the mask
    contents = numpy.vstack([cat[:, RA_col], cat[:, DEC_col],
                            cat[:, z_col], e1_exp, e2_exp, de_exp])
    contents = contents.transpose()
    writettypeCat(contents, headernames, filename)

    return
makeMockCat.__doc__ = \
    """
    Purpose: This is supposed to be a diagnostic test to see how the
    ellipticities should look like given a particular pair of
    M_200 values at the location where the real galaxies are.

    input:
    catalog = string containing filename of the catalog file
                that you want to replicate with a theoretical profile
    M_200_pred = a list of the numerical values of the predicted pair of
                the NFW masses in units of 1e14 M_sun
    coord = a list of strings that denotes what the column names of the
            catalog should be, follows the order e.g. ('ra', 'dec', 'z')
            this follows Will's convention in nfwMCMC.py
    pix_coord = a list of strings that denotes what the column names
                of the catalog should be for pixel coordinates
                e.g. ('x','y')
    ellip_theo = a list of strings that denotes what the column names of
                the theoretical ellipticities predicted by the NFW
                profile and the corresponding uncertainities
                should be in the order of
                e.g. ('e1_theo', 'e2_theo', 'de')
    ellip_inststd = shape noise of the galaxies
    halos = list of numerical numbers, in the form of
            (ra1, dec1, z1, ra2, dec2, z2), RA and DEC positions
            denote the centers of the halos
    r_bounds = a list containing two numbers denoting the inner and
                outer bound of the halo radius.... [in Mpc]
                if it is not specified then it 'll be specified so that
                there is no effects on the halos
    cd = a 2D array that denotes the CD matrix

    output:
    catalog of the input galaxies with the predicted ellipticities
    in pixel coordinates, format follows Will's catalogs
    """

def writettypeCat(contents, headernames, filename, verbose=True):
    '''
    This function
    writes out a ttype catalog compliant with most catalog files used
    by Will Dawson

    Input:
    =====
    contents = numpy array that contains data that you want to write
                out each column of this multi-dimensional numpy array
                should should be one data field.  Note that only numeric
                data can be outputted
    headernames = a list of strings of the names of the variables that you
                 want to write out
    filename = a string that denotes the filename of the catalog

    verbose = logical
        if messages about what the function is doing should be printed out

    Returns:
    =======
    nothing

    One file written out

    Example usage with pandas dataframe:
    $mockCat.writettypeCat(contents = np.array(dataframe),
                            headernames = dataframe.columns,
                            filename = 'filename.txt')

    '''

    # convert the string of headers into ttype headers
    ttypeheaders = process_header(headernames)

    # use the default numpy function in order to write out the catalog
    np.savetxt(filename, contents, delimiter=' ', header=ttypeheaders,
               comments='')
    if verbose == 1:
        print 'ttype catalog named ' + filename + ' has been written'

    return


def addNoise(ellip, sigma_noise):
    '''
    this function adds Gaussian noise with zero mean and sigma =
    sigma_noise to the ellipiticies
    Input:
    ellip = numpy array that contains the data that you want to contaminate

    output:
    ellip with noise
    '''
    noise = [np.random.randn() * sigma_noise for _ in xrange(ellip.size)]

    return ellip + noise



def process_header(headernames):
    ''' read in a list of strings then concatenate the names into the
    ttype header format

    Parameters
    =========
    headernames = a list of strings of the headernames
    if you are using pandas dataframe, just pass in dataframe.columns

    Returns:
    =======
    string - denotes headers formatted correctly for printing / writing out
    for a ttype catalog

    Status:
    ======
    stable and debugged
    '''
    headers = ''
    for i in range(len(headernames)):
        if i == 0:
            headers = '#ttype' + str(i) + ' = ' + headernames[i] + '\n'
        elif i != len(headernames) - 1:
            headers = headers + '#ttype' +\
                str(i) + ' = ' + headernames[i] + '\n'
        elif i == len(headernames) - 1:
            # print 'handling the last header'
            headers = headers + '#ttype' + str(i) + ' = ' + headernames[i]
    return headers


def read_header(file, verbose=True):
    '''
    This functions helps a pandas read_table / read_csv file to read ttype
    headers
    input:
    =====
    file = string,
        contains full path to textfile to be read
    verbose = bool,
        indicates if message should be printed out

    outputs:
    =======
    names = list of strings,
        names of columns that you can use in pandas
            dataframe
    skiprows = integer,
        number of rows that have been skipped

    Stability:
    works, not the fastest nor memory efficient implementation,
    will not work well for large files,
    should probably add exception handling
    should probably use REGEX instead hahaha
    '''
    import re
    f = open(file, 'r')
    l = f.readlines()
    names = []
    skiprows = 0
    for i in range(len(l)):
        a = l[i]
        # strip all whitespace
        a = a.strip()
        # try to only read in the ttype headers
        # this should be written with regex instead in next version
        if a.find('#ttype') == 0 or a.find('# ttype') == 0:
            h = re.split('=', a)
            h = h[1].rstrip('\n')
            h = h.strip()
            # remove all the comments after the header name
            # assuming its a catalog from SEXtractor
            # too adhoc... haha
            k = h.split(' ')
            h = k[0]
            names.append(h)
        # try to find all the lines that are comments
        if a.find('#') == 0:
            skiprows += 1
    if verbose:
        print 'mockCat.read_header:'
        print '# of header parsed = {0}'.format(len(names))
        print 'column names are {0}'.format(names)
        print 'rows to be skipped = {0}'.format(skiprows)

    return names, skiprows


def prepare_fits(fits_header):
    """ get the WCS info from a fits file
    Parameters:
    ===========
    fits_header = string
        denotes the full filepath to the fits_header

    Returns:
    =======
    wcs = astropy wcs object wcs.WCS() containing the header
        see documentation of astropy at
        http://docs.astropy.org/en/latest/wcs/index.html
    cd = a size 2-by-2 list
        ordered like a matrix, [[CD0_0, CD0_1], [CD1_0, CD1_1]]
        and is symmetric in the simpliest cases

    Stability: Works

    """
    hdulist = astropy.io.fits.open(fits_header)
    prihdr = hdulist[0].header
    w = wcs.WCS(prihdr)
    cd = w.wcs.cd
    return w, cd


def make_wcs(cd_matrix):
    """make a wcs object and return it given a cd matrix
    cd_matrix = numpy array of size (2, 2)
    """
    assert cd_matrix.shape == (2, 2), "cd matrix is off wrong shape"
    w = wcs.WCS(naxis=2)
    w.wcs.cd = cd_matrix
    return w


def prepare_halo(halo_pix, w):
    """get the format of pixel values of halo center into wcs form
    Parameters:
    =========
    halo_pix = list of float length in multiples of 3
        with a structure of the 3 entries for each halo to follow
        [x_halo1, y_halo1, z_halo1, xhalo2, y_halo2, z_halo2, ...]
    w = astropy.wcs.wcs.WCS() info

    Returns
    ======
    halo_wcs = list of floats
        with the same dimensions as halo_pix but replaces the pix
        coordinates with RA, DEC then redshift
    halo_num = integer
    """
    assert len(halo_pix) % 3 == 0, "halo_pix list not in multiples of 3"
    assert len(halo_pix) > 0, "halo_pix cannot be empty"
    halo_num = len(halo_pix) // 3

    halos_wcs = []
    for i in range(halo_num):
        # print halo_pix
        x = halo_pix[3 * i]
        y = halo_pix[3 * i + 1]
        z = halo_pix[3 * i + 2]
        [hpix, junk] = w.wcs_pix2world(np.array([[x, y], [x, y]]), 1)
        halos_wcs.append(hpix[0])
        halos_wcs.append(hpix[1])
        halos_wcs.append(z)

    return halos_wcs, halo_num


def prepare_catalog(w, pix_cat, wcs_cat_fullpath, halo_num, e_scale=1.0,
                    bg_gal_z=None, sep=None,
                    pix_col=['x', 'y', 'e1', 'e2', 'de', 'z'], verbose=True):
    """turns the pix coordinates inside a catalog into wcs coordinates
    and outputs the converted catalog to a file

    Parameters
    ==========
    w = astropy.wcs.wcs.WCS object
    pix_cat = string
        full filepath to the ttype catalog with pixel coordinates
    halo_num = float
        number of NFW halos to fit
    bg_gal_z = float
        denotes the background galaxy effective redshift
    pix_col = list of strings
        labels the column header of the output catalog that corresponds to
        the x pixel, y pixel, e1, e2, ellipticity error, redshift

    Returns
    =======
    coord = list of strings for ttype names for the wcs coordinates
    ellip_meas = list of strings for the ellipticities
    """
    # look at the catalog and extract relevant info
    names, skiprows = read_header(pix_cat)

    # read in catalog with a pandas dataframe
    if sep is not None:
        df = \
            pd.read_table(pix_cat, names=names, skiprows=skiprows, sep=sep)
    else:
        df = \
            pd.read_table(pix_cat, names=names, skiprows=skiprows, sep=r'\s*')

    # handle the naming convention of the files
    coord = ['ra', 'dec', 'z']
    ellip_meas = []
    for i in range(1, 3):
        ellip_meas.append('e' + str(i) + '_pix')
        ellip_meas.append('de')

    # use the naming convention to write out the wcs catalog
    if bg_gal_z is not None:
        df[coord[2]] = bg_gal_z * np.ones(df.shape[0])

    # do not forget to scale the ellipticities appropriately
    df[pix_col[2]] = e_scale * df[pix_col[2]].copy(deep=True)
    df[pix_col[3]] = e_scale * df[pix_col[3]].copy(deep=True)


    # rename some of the columns, adhoc solution for this particular
    # catalog, want to make it explicit that we are using ellip. for pixels
    print "WARNING: renaming e1, e2 to e1_pix and e2_pix"
    print "remove line in prepare_catalog if this is not wanted"
    df.rename(columns={pix_col[2]: 'e1_pix', pix_col[3]: 'e2_pix'},
              inplace=True)

    wcs_coord_info = w.wcs_pix2world(df[['x', 'y']], 1)
    df[coord[0]] = wcs_coord_info.transpose()[0]
    df[coord[1]] = wcs_coord_info.transpose()[1]

    # call function to write out the dataframe to file
    writettypeCat(np.array(df), df.columns, wcs_cat_fullpath, verbose)

    return coord, ellip_meas


def plot_ellip(cat, names, ellip_wcs_to_pix, cd=None, sep=r'\s*', scale=1.0,
               addPts_x=[], addPts_y=[], save=False, filename=None,
               showPlot=True, title=None, flipX=False):
    """plot the ellipticites in pixel coordinates
    Parameters
    ==========
    cat = string
        full file path to a galaxy catalog with data specified
        assumes ttype catalog
    sep = string
        denotes the string for separation different columns
        for more details see sep parameter in pandas.read_csv
    names = list of strings
        specifies the names for the headers in the dataframe
        in the order of ['coord1', 'coord2', 'e1_coord', 'e2_coord']
    scale = float
        denotes the scale to multiply the ellipticities by for plotting
        purpose
    ellip_pix_to_wcs = logical
        specifies if the ellipticities in the catalog needs to be
        transformed from pix coordinate to wcs coordinates
    addPts = list
        should have shape = (*, 2)
    save = logical
        specifies if the resulting plot should be saved
    filename = string
        denotes the filepath (name) that should be used for saving
        the ellipticity plot
    verbose = logical
        determines if any messages should be printed

    Returns
    =======
    Nothing?


    Note: by default this assumes that the catalog should be in ttype space
    separated format, if not, modify this function to do what you want
    """
    import cosmoplot as cplot
    import ellip

    colnames, skiprows = read_header(cat)

    if ellip_wcs_to_pix and cd is not None:
        df = pd.read_table(cat, sep=sep, names=colnames, skiprows=skiprows)
        e1_wcs, e2_wcs =\
            ellip.ellip_pix_to_wcs(cd, df[names[2]], df[names[3]])
        names[2] = 'e1_wcs'
        names[3] = 'e2_wcs'
        df[names[2]] = e1_wcs
        df[names[3]] = e2_wcs
    elif ellip_wcs_to_pix and cd is None:
        raise ValueError("Requires a cd matrix")
    else:
        df = pd.read_table(cat, sep=sep, names=colnames, skiprows=skiprows)

    header_exist(df.columns, names)

    if flipX is True:
        cplot.stickplot(-df[names[0]], df[names[1]], df[names[2]],
                        df[names[3]], scale=scale)
    else:
        cplot.stickplot(df[names[0]], df[names[1]], df[names[2]],
                        df[names[3]], scale=scale)

    if ellip_wcs_to_pix:
        xlim_low, xlim_high = plt.xlim()
        plt.xlim(xlim_high, xlim_low)
    plt.ylabel('ellip x {0}'.format(scale), size=15)
    plt.xlabel('ellip x {0}'.format(scale), size=15)
    if title is not None:
        plt.title(title)

    if flipX is True:
        addPts_correct_dimension_and_plot(-np.array(addPts_x), addPts_y)
    else:
        addPts_correct_dimension_and_plot(addPts_x, addPts_y)

    if save and filename is None:
        filename = raw_input("please provide ellip. plot output filename:\n")
        plt.savefig(filename, bbox_inches='tight')
    if showPlot:
        plt.show()

    return


def addPts_correct_dimension_and_plot(addPts_x, addPts_y,
                                           color='red'):
    """ test if list of pts to be plotted has correct dimension then plot
    Paremeters:
    addPts_x = list of float
        x coord of points
    addPts_y = list of float
        y coord of points
    color = string
        color for the points
    """
    if len(addPts_x) + len(addPts_y) > 0:
        assert len(addPts_x) == len(addPts_y),\
            "Shape of the list addPts_x and addPts_y don't match"
        plt.plot(addPts_x, addPts_y, 'x', color=color, mew=2)
    return


def header_exist(colnames, names):
    """test if colnames contain names provided using a dictionary

    Parameters:
    =========
    colnames = list of strings
    names = list of strings
    """
    hashTable = {}
    for key in colnames:
        hashTable[key] = 1

    for key in names:
        try:
            hashTable[key] += 1
        except KeyError:
            raise KeyError("required column name {0} ".format(key) +
                           "missing from catalog header in plot_ellip()")
    return


def double_check_inputs(
        halos_wcs, wcs_cat_fullpath, coord, ellip_meas, ellip_intstd, r_bounds,
        parambounds, N_maxlinks, N_chain, N_burn, N_check, propdist,
        prefix_nfwMCMC, seedbounds, bins, halo_names, cd=((-1, 0), (0, 1)),
        h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0, verbose=True):
    """ask user to manually inspect input variables of nfwMCMC.fitM200
    args = identical to those of fitM200
    kwargs = same as above
    """
    # iterate through all the key-value pairs of the local variables
    for k, v in locals().iteritems():
        sep = "---------------------------------------------------------------"
        print sep
        question = "Does {0} = \n".format(k) + \
            "{0} \n".format(v) + \
            "look ok? (type 'n' for no and quit, type anything else for yes)\n"
        ok = raw_input(question)
        print sep

        if ok.strip() == "n":
            raise ValueError("Input not ok, quitting")

    return


def write_inputs(
        halos_wcs, wcs_cat_fullpath, coord, ellip_meas, ellip_intstd, r_bounds,
        parambounds, N_maxlinks, N_chain, N_burn, N_check, propdist,
        prefix_nfwMCMC, seedbounds, bins, halo_names, cd=((-1, 0), (0, 1)),
        h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0, verbose=True):
    R = open(prefix_nfwMCMC + "_results.txt", "a")
    R.write("--------Input values ---------------\n")
    for k, v in locals().iteritems():
        R.write("{0} = {1}\n".format(k, v))
    R.write("--------end input values -----------\n\n")
    return


def convert_pix2physical(x, y, origin_x, origin_y, wcs):
    return phy_x, phy_y


def compute_shared_var():
    return shared_var


#def compute_data_halo(cat, coord, halo_prop, cosmology, r_bounds, halo_name):
#    metro.data_halo()
#    return data_halo


def compute_ang_sep(cat, coord, halo_prop):
    """
    Parameters
    ==========
    cat = pandas df
        contains the columns named coord
    coord = list of strings
        the names of the columns in the dataframe that corresponds to the
        coords in wcs and in degrees
    halos = list of numbers
        in the form of [RA_wcs_deg, DEC_wcs_deg, redshift, ...]

    Returns
    =======
    numpy array of floats denotes angular separation in units of radians
    """
    return ang_util.angular_separation(cat[coord[0]] / 180. * np.pi,
                                       cat[coord[1]] / 180. * np.pi,
                                       halo_prop._RA / 180. * np.pi,
                                       halo_prop._DEC / 180. * np.pi)


def compute_model_data():
    return


#def write_inputs(*args, **kwargs):
#    halos_wcs, wcs_cat_fullpath, coord, ellip_meas,
#                 ellip_intstd, r_bounds, parambounds, N_maxlinks, N_chain,
#                 N_burn, N_check, propdist, prefix_nfwMCMC, seedbounds,
#                 bins, halo_names, cd=cd, h_scale=h_scale, Om=Om0, Ol=Ol0,
#                 Or=Or0, verbose=verbose):
#
#    for k in args.keys()
#        if k = "prefix_nfwMCMC"
#        prefix = k
#    import pickle
#    f = open("test.pickle")
#    pickle.dump()
#    return
