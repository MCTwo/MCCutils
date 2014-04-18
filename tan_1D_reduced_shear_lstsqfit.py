# coding: utf-8
import numpy as np

# homebrewed modules
from __future__ import division
from nfwMCMC import shear



def tan_1D_reduced_shear(theta, conc, r_s=np.nan):
    '''
    for calculating the tangential shear given by a NFW profile
    for fitting the concentration parameter and the scale radius r_s

    This is written by adopting the expression in Umetsu
    Note that unlike Will 's code this does not depend on the galaxy
    catalog

    Inputs:
    ------
    theta = range of radius of the halo that we 're considering, unit in arcsec
    conc = concentration parameter at the given radius

    r_s = scale radius of the NFW profile (Mpc)
    theta = radius of the halo that is under consideration (arcsec)
    z_halo = halo redshift
    h_scale = hubble scale H = h*100 km/s / Mpc
    Om = matter energy density
    Ol = dark energy density
    Or radiation energy density
    halo coord 1
    halo coord 2
    '''
    # parameters for El Gordo
    z_halo = 0.89
    z_source = 1.22
    beta = 0.216

    # Cosmological parameters
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
    if r_s == np.nan:
        A_200 = 5.71
        B_200 = -0.084
        C_200 = -0.47
        # the h_scale is multiplied because the pivotal mass
        # m_pivotal = 2e12 is in unit of M_sun h_scale^(-1)
        m_200 = (2e12 * h_scale) * \
            (conc / A_200 / (1 + z_halo) ** C_200) * (1 / B_200)
        r_200 = (m_200 / (4 * np.pi / 3 * 200 * rho_cr)) ** (1 / 3.)
        r_s = r_200 / conc

    kappa_s = 2 * del_c * rho_cr * (r_s * minMpc) / nfwSigmacr
    # just do a simple conversion from arcsec to physical distance using
    # angular diameter distance
    # make sure that the units is in arcsec
    theta_s = r_s / dl * 180.0 / np.pi * 60. * 60.
    x = theta / theta_s
    # plt.plot(x)
    # write an array with finer bins than x
    fx = x  # np.arange(np.min(x),np.max(x),(x[1]-x[0])/15.)
    # plt.plot(fx)

    # code the expression for kappa as a function of x
    # print kappa_s
    func_kappa0 = \
        lambda fx: kappa_s / \
        (1 - fx ** 2.) * (-1. + 2. / np.sqrt(1. - fx ** 2.)
                          * np.arctanh(np.sqrt(1. - fx) / np.sqrt(1. + fx)))
    func_kappa1 = lambda fx: kappa_s * 1. / 3.
    func_kappa2 = lambda fx: kappa_s / \
        (fx ** 2. - 1) * \
        (1. - 2. / np.sqrt(fx ** 2. - 1.)
         * np.arctan(np.sqrt(fx - 1.) / np.sqrt(fx + 1.)))

    kappa = np.piecewise(fx, [fx < 1.0, fx == 1.0, fx > 1.0],
                         [func_kappa0, func_kappa1, func_kappa2])

    g_theta = \
        np.piecewise(x, [x < 1.0, x == 1.0, x > 1.0],
                     [lambda x: np.log(x / 2.) + 2. / np.sqrt(1 - x ** 2)*\
                      np.arctanh(np.sqrt(1 - x) / np.sqrt(1 + x)),
                      lambda x: np.log(x / 2.) + 1.,
                      lambda x: np.log(x / 2.) + 2. /\
                      np.sqrt(x ** 2. - 1.) * np.arctan(np.sqrt(x - 1.) /\
                                                        np.sqrt(x + 1.))])
    kappa_bar = 2. * kappa_s / x ** 2. * g_theta
    # plt.plot(fx*theta_s,kappa)
    # we did not have any angular info for kappa to begin with
    azim_kappa = kappa
    #azim_bins = azimuthal_avg_shear_in_bins(kappa, fx, np.min(x), np.max(x)+x[1]-x[0], x[1]-x[0])

    red_shear = (kappa_bar - azim_kappa) / (1 - azim_kappa)

    return red_shear
