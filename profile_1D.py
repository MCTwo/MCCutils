#-------------------------------------------------------------------------
# Author: Karen Yin-Yee Ng <karenyng@ucdavis.edu>
# Purpose for plotting 1D reduced shear profile of a NFW halo
# Date: 07/04/2013
# License: BSD
#-------------------------------------------------------------------------

from __future__ import division
import cosmo
import profiles
import numpy as np


def azimuthal_avg_shear_in_bins(k, x, startbin=0, endbin=10, b_width=1):
    '''
    Purpose:
    average the value of shear within a certain radius
    Input:
    k = shear that you want to bin
    startbin = the value at the starting bin
    endbin = the value at the ending bin
    b_width = the width of the bin
    '''
    bins = int((endbin - startbin) / b_width)
    # print '# of bins is ',bins
    # initialize bin edge
    bin_edge = np.arange(startbin, endbin + b_width, b_width)
    # print 'bin_edge is ',bin_edge
    # initialize the average value of each bin
    bin_array = np.arange(
        startbin +
        b_width /
        2.0,
        endbin +
        b_width /
        2,
        b_width)
    # print 'bin array is', bin_array
    # exclude values outside the desired range
    mask = np.logical_and(x > startbin, x <= endbin)
    red_x = x[mask]
    red_k = k[mask]

    azim_kappa = np.ones(bin_array.size)

    for i in range(bins):
        # only look at values within each bin within each iteration
        temp_mask = np.logical_and(
            red_x >= bin_edge[i],
            red_x < bin_edge[i + 1])
        if red_k[temp_mask].size == 0:
            temp_k = 0
        else:
            temp_k = red_k[temp_mask]
        # average the e1_pix value within this bin
        if np.isnan(np.mean(temp_k)):
            raise ValueError('NaN encountered')
        azim_kappa[i] = np.mean(temp_k)
        # print '# of data point in this bin = ',len(temp_k), ',bin =
        # ',(bin_edge[i]+bin_edge[i+1] )/2., ' bin_array=',bin_array[i]
    return azim_kappa, bin_array


def ext_1D_reduced_shear(theta, conc, r_s):
    '''
    function for calculating the tangential shear given by a NFW profile
    for fitting the concentration parameter and the scale radius r_s

    This is written by adopting the expression in Umetsu

    theta = range of radius of the halo that we 're considering,
    unit in arcsec
    conc = concentration parameter at the given radius

    r_s = scale radius of the NFW profile (Mpc)
    z_halo = halo redshift

    ### WARNING
    ### cosmological parameters are written within the function

    h_scale = hubble scale H = h*100 km/s / Mpc
    Om = matter energy density
    Ol = dark energy density
    Or = radiation energy density
    halo coord 1
    halo coord 2

    z_halo = redshift of the lens
    z_source = redshift of the source galaxies
    beta = ratio of D_LS and D_S see James Jee 's paper for exact def

    '''
    # latest parameters for El Gordo
    z_halo = 0.87
    z_source = 1.318
    beta = 0.258

    # old parameters
    #z_halo =0.89
    #z_source = 1.22
    #beta = 0.216

    # Cosmological parameters
    # print "profile_1D.tan_1D_reduced_shear: this func uses its own "
    # print "cosmological parameters"
    Om = 0.3
    Ol = 0.7
    Or = 0.0
    c = 3e5  # speed of light units km/s
    G = 6.673 * 10 ** (-11)  # m^3/kg/s^2
    h_scale = 0.7
    kminMpc = 3.08568025 * 10 ** 19  # km in a Megaparsec
    minMpc = 3.08568025 * 10 ** 22  # m in a Megaparsec
    kginMsun = 1.988e30  # kg

    # angular diameter distance to the lens
    dl = cosmo.Da(z_halo, h_scale, Om, Ol)  # in Mpc
    del_c = 200 / 3. * conc ** 3 / \
        (np.log(1 + conc) - conc / (1 + conc))  # unitless
    # in units of kg/m^2
    nfwSigmacr = c ** 2 / (4 * np.pi * G * dl * beta) * 1000 / kminMpc
    # print 'nfwSigmacr in units of M_sun /pc^2 is ',
    # nfwSigmacr/kginMsun*(3.086*10**16)**2

    # in units of kg/m^3
    rho_cr = cosmo.rhoCrit(z_halo, h_scale, Om, Ol, Or)  # /kginMsun*minMpc**3

    # if scale radius is not given then infer from Duffy et al.
    # assuming full halo

    kappa_s = 2 * del_c * rho_cr * (r_s * minMpc) / nfwSigmacr
    # just do a simple conversion from arcsec to physical distance using
    # angular diameter distance
    # make sure that the units is in arcsec
    theta_s = r_s / dl * 180.0 / np.pi * 60. * 60.
    x = theta / theta_s

    # write an array with finer bins than x
    fx = x  # np.arange(np.min(x),np.max(x),(x[1]-x[0])/10.)

    # code the expression for kappa as a function of x
    # print kappa_s
    func_kappa0 = lambda fx: kappa_s / (1 - fx ** 2.) *\
        (-1. + 2. / np.sqrt(1. - fx ** 2.) *
         np.arctanh(np.sqrt(1. - fx) / np.sqrt(1. + fx)))
    func_kappa1 = lambda fx: kappa_s * 1. / 3.
    func_kappa2 = lambda fx: kappa_s / (fx ** 2. - 1) *\
        (1. - 2. / np.sqrt(fx ** 2. - 1.) *
         np.arctan(np.sqrt(fx - 1.) / np.sqrt(fx + 1.)))

    kappa = np.piecewise(fx, [fx < 1.0, fx == 1.0, fx > 1.0],
                         [func_kappa0, func_kappa1, func_kappa2])

    g_theta = np.piecewise(x, [x < 1.0, x == 1.0, x > 1.0],
                           [lambda x: np.log(x / 2.) + 2. / np.sqrt(1 - x ** 2) *
                            np.arctanh(np.sqrt(1 - x) / np.sqrt(1 + x)),
                            lambda x: np.log(x / 2.) + 1.,
                            lambda x: np.log(x / 2.) + 2. / np.sqrt(x ** 2. - 1.) *
                            np.arctan(np.sqrt(x - 1.) / np.sqrt(x + 1.))])
    kappa_bar = 2. * kappa_s / x ** 2. * g_theta
    azim_kappa = kappa
    # no need to azimuthally smooth function since
    # the function did n't contain azimuthal angular info to begin with
    #, azim_bins = azimuthal_avg_shear_in_bins(kappa, fx, np.min(x), np.max(x)+x[1]-x[0], x[1]-x[0])

    red_shear = (kappa_bar - azim_kappa) / (1 - azim_kappa)

    return red_shear


def tan_1D_reduced_shear(theta, z_halo, z_source, beta, conc):
    '''
    function for calculating the tangential shear given by a NFW profile
    for fitting the concentration parameter and the scale radius r_s

    This is written by adopting the expression in Umetsu

    theta = numpy vector of floats
        range of radius of the halo that we 're considering, unit in arcsec
    z_halo = float, halo redshift
    z_source = float, effective source galaxy redshifts
    beta = float, ratio of D_LS and D_S see James Jee 's paper for exact def
    conc = float, concentration parameter at the given radius


    ### WARNING
    ### cosmological parameters are written within the function

    h_scale = hubble scale H = h*100 km/s / Mpc
    Om = matter energy density
    Ol = dark energy density
    Or radiation energy density
    '''
    # old parameters
    #z_halo =0.89
    #z_source = 1.22
    #beta = 0.216

    # Cosmological parameters
    # print "profile_1D.tan_1D_reduced_shear: this func uses its own "
    # print "cosmological parameters"
    Om = 0.3
    Ol = 0.7
    Or = 0.0
    c = 3e5  # speed of light units km/s
    G = 6.673 * 10 ** (-11)  # m^3/kg/s^2
    h_scale = 0.7
    kminMpc = 3.08568025 * 10 ** 19  # km in a Megaparsec
    minMpc = 3.08568025 * 10 ** 22  # m in a Megaparsec
    kginMsun = 1.988e30  # kg

    # angular diameter distance to the lens
    dl = cosmo.Da(z_halo, h_scale, Om, Ol)  # in Mpc
    del_c = 200 / 3. * conc ** 3 / \
        (np.log(1 + conc) - conc / (1 + conc))  # unitless
    # in units of kg/m^2
    nfwSigmacr = c ** 2 / (4 * np.pi * G * dl * beta) * 1000 / kminMpc
    # print 'nfwSigmacr in units of M_sun /pc^2 is ',
    # nfwSigmacr/kginMsun*(3.086*10**16)**2

    # in units of kg/m^3
    rho_cr = cosmo.rhoCrit(z_halo, h_scale, Om, Ol, Or)  # /kginMsun*minMpc**3

    # if scale radius is not given then infer from Duffy et al.
    # assuming full halo
    # if r_s == np.nan:
    r_s = scale_radius(conc, z_halo, Om, Ol, Or)  # in Mpc

    kappa_s = 2 * del_c * rho_cr * (r_s * minMpc) / nfwSigmacr
    # just do a simple conversion from arcsec to physical distance using
    # angular diameter distance
    # make sure that the units is in arcsec
    theta_s = r_s / dl * 180.0 / np.pi * 60. * 60.
    x = theta / theta_s

    # write an array with finer bins than x
    fx = x  # np.arange(np.min(x),np.max(x),(x[1]-x[0])/10.)

    # code the expression for kappa as a function of x
    # print kappa_s
    func_kappa0 = lambda fx: kappa_s / (1 - fx ** 2.) *\
        (-1. + 2. / np.sqrt(1. - fx ** 2.) *
         np.arctanh(np.sqrt(1. - fx) / np.sqrt(1. + fx)))
    func_kappa1 = lambda fx: kappa_s * 1. / 3.
    func_kappa2 = lambda fx: kappa_s / (fx ** 2. - 1) *\
        (1. - 2. / np.sqrt(fx ** 2. - 1.) *
         np.arctan(np.sqrt(fx - 1.) / np.sqrt(fx + 1.)))

    kappa = np.piecewise(fx, [fx < 1.0, fx == 1.0, fx > 1.0],
                         [func_kappa0, func_kappa1, func_kappa2])

    g_theta = \
        np.piecewise(x, [x < 1.0, x == 1.0, x > 1.0],
                     [lambda x: np.log(x / 2.) + 2. / np.sqrt(1 - x ** 2) *
                      np.arctanh(np.sqrt(1 - x) / np.sqrt(1 + x)),
                      lambda x: np.log(x / 2.) + 1.,
                      lambda x: np.log(x / 2.) + 2. / np.sqrt(x ** 2. - 1.) *
                      np.arctan(np.sqrt(x - 1.) / np.sqrt(x + 1.))])

    kappa_bar = 2. * kappa_s / x ** 2. * g_theta
    azim_kappa = kappa  # no need to azimuthally smooth function since
    # the function did n't contain azimuthal angular info to begin with
    #, azim_bins = azimuthal_avg_shear_in_bins(kappa, fx, np.min(x), np.max(x)+x[1]-x[0], x[1]-x[0])

    # the units do not look right
    red_shear = (kappa_bar - azim_kappa) / (1 - azim_kappa)

    return red_shear


def scale_radius(conc, z_halo, Om=0.3, Ol=0.7, Or=0.0,
                 A_200=5.71, B_200=-0.084, C_200=-0.47, h_scale=0.7):
    '''
    purpose: compute the scale radius given the concentration parameter
    default parameters based on full profiles of Duffy et al. 2008
    input:
    conc = concentration parameter at r200
    z_halo = redshift of that halo
    Om, Ol, Or = cosmological paramters
    A200 = duffy et al parameter for concentration radius relationship
    B200 = "
    C200 = "
    h_scale = reduced hubble parameter

    output:
    scale radius in Mpc
    '''
    # unit conversion values
    minMpc = 3.08568025 * 10 ** 22  # m in a Megaparsec
    kginMsun = 1.988e30  # kg

    rho_cr = cosmo.rhoCrit(z_halo, h_scale, Om, Ol, Or)  # in kg/m**3
    # the h_scale is multiplied because the pivotal mass
    # m_pivotal = 2e12 is in unit of M_sun h_scale^(-1)
    m_200 = profiles.nfwM200(conc, z_halo, A_200, B_200, C_200, h_scale) *\
        kginMsun
    r_200 = (m_200 / (4 * np.pi / 3 * 200 * rho_cr)) ** (1 / 3.)  # in m
    r_s = r_200 / conc / minMpc

    return r_s


class halo:
    """
    group halo properties together for fitting

    protected variables:
    -------------------
    self._conc = float, concentration of the NFW halo
    self._R200c = float, in unit of Mpc
    self._M200c = float, in units of 1e14 Msun
    self._comments = comment about what the values of this halo is from
    """
    def __init__(self, conc, R200c, M200c, comments=None):
        self._conc = conc
        self._R200c = R200c
        self._M200c = M200c
        self._comments = comments


