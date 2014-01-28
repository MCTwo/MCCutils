'''
Author: Karen Y. Ng (karenyng@ucdavis.edu)
Date: 10/05/2013
License: 3 clause BSD
Note:
To be used in conjunction with the functions in nfwMCMC.py by Will Dawson
'''
from __future__ import division
#import pdb
import numpy
import numpy as np
#import pylab
import cosmo
import tools
#from profiles import nfw_Sigma
#from profiles import nfw_Sigmabar
from profiles import nfwparam
import numpy.random
from nfwMCMC import shear
from astropy import wcs

def makeMockCat(catalog, M_200_pred, coord, ellip_theo, ellip_inststd,
                halos, filename, r_bounds=(0, 1000), noise_sigma=0.0,
                pix_coord=('x','y'), h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0,
                cd = ((-1,0),(0,1)), verbose=False):
    '''
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
    '''

    #read in the catalogs
    cat = tools.readcatalog(catalog)
    key = tools.readheader(catalog)
    #number of galaxies
    N_gal = numpy.shape(cat)[0]

    #-----------Perform some basic user input checks----------------
    #Determine the number of halo's to model
    N_halo = numpy.size(halos)//3
    if numpy.size(halos)%3. != 0:
        print 'fitM200: an invalid "halos" input was specified, exiting'
        sys.exit()

    if verbose:
        print 'fitM200: calculating properties of the galaxies'+\
            ' with respect to each halo'

    #-----------Initialize some of the vectors ---------------------
    del_a = numpy.zeros((N_gal,N_halo))
    del_d = numpy.zeros((N_gal,N_halo))
    invSigmacr = numpy.zeros((N_gal,N_halo))
    #var_invSigmacr = numpy.zeros((N_gal,N_halo))
    mask_r_inner = numpy.zeros((N_gal,N_halo)) != 0
    mask_r_outer = numpy.zeros((N_gal,N_halo)) != 0
    mask_z = numpy.zeros((N_gal,N_halo)) != 0

    for h in numpy.arange(N_halo):
        if verbose:
            print 'fitM200: evaluating galaxies with respect'+\
                    ' to halo {0}'.format(h)
        for g in numpy.arange(N_gal):
            #Calculate the angular separations
            del_a[g,h],del_d[g,h] = tools.angcomp(cat[g,key[coord[0]]],
                                                  cat[g,key[coord[1]]],
                                                  halos[h*3],halos[h*3+1])
            #*# Note that I could speed up the code by masking galaxies
            #outside
            #*# the apodizing radii and not calculating their invSigmacr
            #Calulate inverse sigma critical
            if halos[h*3+2] < cat[g,key[coord[2]]]:
    #            print 'right before invsigcr calc'
    #            print 'z_halo = {0}'.format(halos[h*3+2])
    #            print 'z_gal = {0}'.format(cat[g,key[coord[2]]])
                invSigmacr[g,h] = 1/cosmo.lensgeo(halos[h*3+2],
                                                  cat[g,key[coord[2]]],
                                                  h_scale,Om,Ol)['sigcr']
    #            print 'right after invsigcr calc'
    #            print 'break point'
            #Calculate the standard error on sigma critcal
            ###################################################
            ### Need to add this section
            ###################################################
        #convert deg anglular separation to Mpc separation
        del_a[:,h] *= 60*cosmo.ProjectedLength(halos[h*3+2],h_scale,Om,Ol)
        del_d[:,h] *= 60*cosmo.ProjectedLength(halos[h*3+2],h_scale,Om,Ol)
    del_r = numpy.sqrt(del_a**2+del_d**2)
    for h in numpy.arange(N_halo):
        #Mask galaxies that are outside the apodizing radii or
        #foreground of the halo
        mask_r_inner[:,h] = del_r[:,h] > r_bounds[0]
        #mask out galaxies that are not within any outer radii of a halo
        mask_r_outer[:,h] = del_r[:,h] < r_bounds[1]
    #    mask_r[:,h] = numpy.logical_and(del_r[:,h]>r_bounds[0],
    #                                    del_r[:,h]<r_bounds[1])
        #this masks the galaxies that are in front of the
        mask_z[:,h] = cat[:,key[coord[2]]] > halos[h*3+2]
    #mask out all the galaxies that are within any of the inner radii
    #of either halo
    mask_r_inner = numpy.sum(mask_r_inner,axis=1) == N_halo

    #mask = mask_r*mask_z
    mask = numpy.reshape(mask_r_inner,(numpy.size(mask_r_inner),1))*\
            mask_r_outer*mask_z

    #mask = mask_r*mask_z
    mask = numpy.reshape(mask_r_inner,(numpy.size(mask_r_inner),1))*\
           mask_r_outer*mask_z

    if verbose:
        for h in numpy.arange(N_halo):
            print 'fitM200: There are {0}'.format(numpy.sum(mask[:,h]))+\
                ' background galaxies for Halo {1}'.format(h)

    ############ there are some initialization and calculation missing
    ####### here !!!!!!

    #Find the unit relation between del_x, del_y and del_a, del_d

    cd_inv = numpy.linalg.inv(numpy.array(cd))
    cd_mag = numpy.sqrt(numpy.sum(cd_inv[0,:]**2))
    del_x = cd_inv[0,0]/cd_mag*del_a+cd_inv[0,1]/cd_mag*del_d
    del_y = cd_inv[1,0]/cd_mag*del_a+cd_inv[1,1]/cd_mag*del_d
    #------Calculate the  angle of the galaxy CCW from +x axis of the halo
    #-------use trig identities to ---------------------
    #-------sinphi = del_y/del_r -----------------------
    #-------cos2phi = 1-2sinphi**2 ---------------------
    #-------sin2phi = 2sinphi*cosphi -------------------
    cos2phi = 1-2*(del_y/del_r)**2
    sin2phi = 2*del_x*del_y/del_r**2

    N_halo = numpy.shape(del_r)[1]
    N_gal = numpy.shape(del_r)[0]
    gamma = numpy.zeros((N_gal,N_halo))
    z_halo = numpy.zeros(N_halo)
    del_c = numpy.zeros(N_halo)
    r_s = numpy.zeros(N_halo)
    for h in numpy.arange(N_halo):
        z_halo[h] = halos[h*3+2]
        #print z_halo[h]

    # Calculate the NFW parameters corresponding to M_200
    #print 'nfwMCMC.nfwparam - M_200 is ',M_200
    for h in numpy.arange(N_halo):
        del_c[h], r_s[h] =nfwparam(M_200_pred[h],z_halo[h],h_scale=0.7,
                            Om=0.3,Ol=0.7,Or=0.0)
        #loop through each halo calculating the expected absolute shear
        #due to that individual halo
        #*# Could place a mapping based multiprocessing function here
        #calculate the expected ellipticity components of each galaxy
        gamma[mask[:,h],h] = shear(del_c[h],r_s[h],del_r[mask[:,h],h],
                             z_halo[h],invSigmacr[mask[:,h],h], h_scale,
                             Om,Ol,Or)
    mask_e = numpy.sum(mask,axis=1) != 0
    #--calculate the expected ellipticities of each galaxy due to all halos
    e1_exp = numpy.sum(-gamma[mask_e,:]*cos2phi[mask_e,:],axis=1)
    e2_exp = numpy.sum(-gamma[mask_e,:]*sin2phi[mask_e,:],axis=1)

    e1_exp = addNoise(e1_exp, noise_sigma)
    e2_exp = addNoise(e2_exp, noise_sigma)


    #need to check the dimensionality of e1_exp and e2_exp after applying
    #the masks

    #---write out the ellipiticies and the corresponding RA / DEC / X / Y
    #---coordinate to a ttype header list
    headernames = (coord[0], coord[1],coord[2], ellip_theo[0],
                    ellip_theo[1], ellip_theo[2])
    #### this is the way to call a column in the catalog e.g.
    #### cat[:,key['name']



    #----------figure out all the key names --------
    RA_col = key[coord[0]]
    DEC_col = key[coord[1]]
    z_col = key[coord[2]]
    #    e1_col = ellip_theo[0]
    #    e2_col = ellip_theo[1]
    #    de_col = ellip_theo[2]

    #-------initialize an array that denotes std dev of the ellipticities
    de_exp = np.zeros(cat.shape[0])

    #stack all the catalog content together in one big 2D numpy array
    #remember to apply the mask
    contents = numpy.vstack([ cat[:,RA_col], cat[:,DEC_col],
                            cat[:,z_col], e1_exp, e2_exp, de_exp ])
    contents = contents.transpose()
    writettypeCat(contents, headernames, filename)

    return


def writettypeCat(contents, headernames, filename, verbose = True):
    '''
    This function
    writes out a ttype catalog compliant with most catalog files used
    by Will Dawson

    Input:
    contents = numpy array that contains data that you want to write
                out each column of this multi-dimensional numpy array
                should should be one data field.  Note that only numeric
                data can be outputted
    headernames = a list of strings of the names of the varaibles that you
                 want to write out
    filename = a string that denotes the filename of the catalog

    Output: return nothing

    One file written out

    Example usage with pandas dataframe:
    $mockCat.writettypeCat(contents = np.array(dataframe),
                            headernames = dataframe.columns,
                            filename = 'filename.txt')

    '''

    #convert the string of headers into ttype headers
    ttypeheaders = process_header(headernames)

    #use the default numpy function in order to write out the catalog
    np.savetxt(filename, contents, delimiter = ' ', header = ttypeheaders,
            comments='')
    if verbose == True:
        print 'ttype catalog named '+ filename + ' has been written'

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
    noise = [np.random.randn()*sigma_noise for _ in xrange(ellip.size)]

    return ellip + noise


def process_header(headernames):
    '''
    Purpose: read in a list of strings then concatenate the names into the ttype header format
    Parameters
    =========
    headernames = a list of strings of the headernames
    if you are using pandas dataframe, just pass in dataframe.columns

    Status:
    ======
    stable and debugged
    '''
    headers = ''
    for i in range(len(headernames)):
        if i == 0:
            headers = '#ttype'+str(i)+' = '+headernames[i]+'\n'
        elif i !=len(headernames)-1:
            headers = headers + '#ttype'+str(i)+' = '+headernames[i]+'\n'
        elif i == len(headernames)-1:
            #print 'handling the last header'
            headers = headers + '#ttype'+str(i)+' = '+headernames[i]
    return headers



def read_header(file, verbose = True):
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
        #try to find all the lines that are comments
        if a.find('#') == 0:
            skiprows += 1
    if verbose:
        print 'mockCat.read_header:'
        print '# of header parsed = {0}'.format(len(names))
        print 'rows to be skipped = {0}'.format(skiprows)

    return names, skiprows


#class cosmo_param:
#    '''
#    Purpose: group all the cosmological parameters together and
#    just pass it as a single object whenever we need to know the cosmology
#    it is too inefficient to pass the parameters one by one
#    '''
#    def Ol(self):
#       '''
#       This represents Omega Lambda
#       '''
#
#    def Om(self):
#
#    def Or(self):


def prepare_fits(fits_header):
    """ get the WCS info from a fits file
    Parameters:
    ===========
    fits_header = string
        denotes the path to the fits_header

    Returns:
    =======
    wcs = astropy wcs object wcs.WCS() containing the header
        see documentation of astropy at
        http://docs.astropy.org/en/latest/wcs/index.html
    cd = a size 2-by-2 list
        ordered like a matrix, [[CD0_0, CD0_1], [CD1_0, CD1_1]]

    Stability: Untested
    =========

    """
    hdulist = astropy.io.fits.open(fits_header)
    prihdr = hdulist[0].header
    w = wcs.WCS(prihdr)

    cd = [[ w.cd[0,0], w.cd[0,1]], [ w.cd[1,0], w.cd[1,1]]
    return w, cd
