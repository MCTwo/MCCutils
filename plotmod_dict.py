"""
This program takes the output pickle arrays of TSM.py and creates some plots
and statistics.  It is largely based on plotTSM.py

The main function of this program is to generate report quality plots,
especially a covariance array plot
"""
from __future__ import division
import pylab
import numpy as np
import numpy
import pickle
from scipy.stats import norm
from scipy.special import erf
import pdb
from astrostats import biweightLoc, bcpcl
import cosmo
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd


def load_pickles_to_df(par, prefix, index, msun1e14=True,
                       msun_string=['m_1', 'm_2'], verbose=False):
    """
    Parameters
    =========
    par = list of strings
        denotes the name of the parameters
    prefix = string
        prefix that denotes the path and file name excluding extension
        .pickle
    index = number
        how many pickle files there are
    msun1e14 = logical
        if we want units of mass to be in 1e14 M_sun
    msun_string = list of strings
        the strings are the keys to the df for changing the units
    verbose: logical

    Returns
    ======
    Pandas Dataframe that have all nans removed

    """
    for i in range(len(par)):
        if(verbose):
            print 'loading ' + par[i]
        d = loadcombo(prefix, index, par[i])
        if i == 0:
            data = pd.DataFrame(d, columns=[par[i]])
        else:
            data[par[i]] = pd.DataFrame(d)

    print "dropping NA!"
    data = data.dropna()

    if msun1e14 is True:
        for mass in msun_string:
            print "converting entry " + mass + \
                " to units of 1e14 m_sun"
            data[mass] = data[mass] / (1e14)

    return data


def loadcombo(prefix, index, suffix):
    """loads the data from pickle files
    Parameters:
    ===========
    prefix = string
        denotes the name of the path to the file
    index = number
        denotes the index to add to the file name
    filename consists of prefix + suffix+.pickle see below

    Returns
    =======
    numpy array
        data
    """
    array = []
    for i in index:
        filename = prefix + i + '_' + suffix + '.pickle'
        # read in the pickled array
        F = open(filename)
        tmp = pickle.load(F)
        F.close()
        array = numpy.append(array, tmp)

    #filename = prefix+suffix+'.pickle'
    # read in the pickled array
    #F = open(filename)
    #tmp = pickle.load(F)
    # F.close()
    #array = numpy.append(array,tmp)
    return array


# Create the prior from the radio relic constraints
# we know that if the two subclusters are within 0.5 Mpc to 1.5 Mpc
# then the detection of a radio relic is possible
# d_3d is in Mpc
# return a mask that can be applied to other data arrays
# also return how many non-zero entries
# def radio_dist_prior(d_3D, mask, d_3Dmax = 3.0, d_3Dmin = 1.0):
#    '''
#    to be modified such that this takes in the range for the uniform prior
#    input:
#    d_3D = numpy array to be masked
#    radiomask = numpy array that is contains either 1 or 0
#    d_3Dmax = float, the upper limit to be masked out
#    d_3Dmin = float, the lower limit to be masked out
#    '''
#    count = 0
#    for r in range(len(d_3D)):
#        if (d_3D[r]>d_3Dmax)  or (d_3D[r]< d_3Dmin):
#            radiomask[r] = 0
#            count += 1
#    count = len(radiomask) - count
#    return radiomask, count
def radio_dist_prior(d_3D, d_3Dmax=3.0, d_3Dmin=1.0):
    '''
    Stability: to be tested
    input:
    d_3D = numpy array to be masked, in unit of Mpc
    d_3Dmax = float, the upper limit to be masked out, in unit of Mpc
    d_3Dmin = float, the lower limit to be masked out, in unit of Mpc
    output:
    mask = numpy array that gives 1 if it is within the range
            0 if it is NOT within the specified range
    count = number of entries along the array that has value 1
    '''
    mask = np.logical_and(d_3D < d_3Dmax, d_3D >= d_3Dmin)
    count = np.sum(mask)
    return mask, count


def radio_polar_prior(alpha, alpha_min=0., alpha_max=40):
    '''
    Stability: to be tested
    input:
    alpha = numpy array to be masked, in units of degrees
    alpha_min = float, if alpha is smaller than this value it's masked out
    alpha_max = float, if alpha is bigger than this value, it's masked out
    output:
    mask = numpy array that gives 1 if it is within the range and 0
            otherwise
    count = number of entries that should remain after masking
    '''
    mask = np.logical_and(alpha < alpha_max, alpha > alpha_min)
    count = np.sum(mask)
    return mask, count


def apply_radioprior(radiomask, dataarray):
    '''
    Checks if the length data array is the same as the prior mask
    if not, do nothing
    if lengths are the same, apply the prior
    starts examine if the length of mask is the same
    as the length of the array to be masked
    input:
    mask = numpy array with true or false as the values
    dataarray = numpy data array to be masked
    '''
    if len(radiomask) != len(dataarray):
        print 'length of mask and data array does not match!'
        print 'skipping the application of radio relic prior'
    else:
        # apply the mask
        temp = dataarray * radiomask

        counter = 0
        # removed the entries that are zero from the data array
        # to avoid the zero entries being binned
        for n in range(len(temp)):
            if temp[n] != 0.0:
                counter += 1
        # print 'number of non-zero entries after masking is', counter

        dataarray = numpy.zeros(counter)
        ncounter = 0
        for n in range(len(temp)):
            if temp[n] != 0.0:
                dataarray[ncounter] = temp[n]
                ncounter += 1
        # print 'number of non-zero entries for data array after masking is',
        # ncounter

    return dataarray

# this function takes in parameter arrays
# bins the data then calculates the percentage difference


def percentdiff(x, prefix, prob=None, N_bins=100, histrange=None, x_lim=None,
                y_lim=None, x_label=None, y_label=None, legend=None):

    #fig = pylab.figure()
    # find out size of array
    totalsize = len(x)
    print 'data size of each variable is ', totalsize
    # divide total number of data points into nparts-1
    nparts = 101
    d = (nparts - 1, 5)
    reduced_d = (nparts-2, 5)
    data = numpy.zeros(d)
    x_perdiff = numpy.zeros(reduced_d)

    # iterate from 1 to nparts-1
    for n in range(1, nparts):
            # print i,"th iteration"
        # size of each of the n parts:
        partsize = totalsize*(n)/(nparts-1)
        hist, binedges, tmp = pylab.hist(
            x[:partsize], bins=N_bins, histtype='step',                              weights=prob[:partsize], range=histrange, color='k',
            linewidth=2)
        # Calculate the location and confidence intervals
        # Since my location and confidence calculations can't take weighted data I
        # need to use the weighted histogram data in the calculations
        for i in numpy.arange(N_bins):
            if i == 0:
                x_binned = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
            else:
                x_temp = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
                x_binned = numpy.concatenate((x_binned, x_temp))
        # print 'len of x_binned is ',len(x_binned)
        loc = biweightLoc(x_binned)
        ll_68, ul_68 = bcpcl(loc, x_binned, 1)
        ll_95, ul_95 = bcpcl(loc, x_binned, 2)
        # print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(prefix,loc,ll_68,ul_68,ll_95,ul_95)
            # this will store the data
        data[n-1] = [loc, ll_68, ul_68, ll_95, ul_95]

        '''
        if n%10 == 0:
            # Create location and confidence interval line plots
            pylab.plot((loc,loc),(pylab.ylim()[0],pylab.ylim()[1]),'--k',linewidth=2,label='$C_{BI}$')
            pylab.plot((ll_68,ll_68),(pylab.ylim()[0],pylab.ylim()[1]),'-.',
            linewidth=2,color='#800000',label='68% $IC_{B_{BI}}$')
            pylab.plot((ul_68,ul_68),(pylab.ylim()[0],pylab.ylim()[1]),'-.',linewidth=2,color='#800000')
            pylab.plot((ll_95,ll_95),(pylab.ylim()[0],pylab.ylim()[1]),':',linewidth=2,color='#0000A0',label='95% $IC_{B_{BI}}$')
            pylab.plot((ul_95,ul_95),(pylab.ylim()[0],pylab.ylim()[1]),':',linewidth=2,color='#0000A0')

            if x_label != None:
                pylab.xlabel(x_label,fontsize=14)
            if y_label != None:
                pylab.ylabel(y_label,fontsize=14)
            if x_lim != None:
                pylab.xlim(0,x_lim)
            if y_lim != None:
                pylab.ylim(0,y_lim)
            if legend != None:
                pylab.legend()
            filename = prefix+'_histplot1D_'+str(n)
            pylab.savefig(filename)
            pylab.close()
        '''
    filename = prefix+'_histplot1D_percentdiff.png'
    pylab.plot((loc, loc), (pylab.ylim()[0], pylab.ylim()[1]),
               '--k', linewidth=2, label='$C_{BI}$')
    pylab.plot((ll_68, ll_68), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#800000', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68, ul_68),
               (pylab.ylim()[0], pylab.ylim()[1]), '-.', linewidth=2, color='#800000')
    pylab.plot((ll_95, ll_95), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#0000A0', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95, ul_95),
               (pylab.ylim()[0], pylab.ylim()[1]), ':', linewidth=2, color='#0000A0')

    pylab.savefig(filename, dpi=300, bbox_inches='tight')
    pylab.close()

    print '\n'+prefix+' data is '
    print data
    print '      '

    for n in range(1, nparts-1):
        x_perdiff[n-1] = (data[nparts-2]-data[n-1])*2 * \
            100/(data[nparts-2]+data[n-1])
    print prefix+' per diff is '
    print x_perdiff
    print '      '

    # this will invert the array, now disabled
    #x_perdiff = x_perdiff[::-1]

    return x_perdiff


# plot the pdf of the histograms
# might want to return the pdf later on
def histplot1d_pdf(x, prefix, prob=None, N_bins=100, histrange=None,
                   x_lim=None, y_lim=None, x_label=None, y_label=None,
                   legend=None, title=None, save=True, verbose=True):
    """plot data as normalized pdf
    x = numpy array like object, can be dataframe columns
        data to be plotted on the x-axis
    prefix = string
        denotes the output file prefix
    prob = numpy array with same size as x
        denotes the weight to be put for correcting bias
    N_bins = integer
        denotes the number of bins
    histrange = a size 2 numpy array / list
        denotes the lower and upper range for making the histogram
    x_lim = a size 2 numpy array
    y_lim = a size 2 numpy array
    x_label = string
    y_label = string
    legend = string
    title = string
    """
    hist, binedges, tmp = \
        pylab.hist(x, bins=N_bins, histtype='step',
                   weights=prob, range=histrange, color='k', linewidth=2)

    # do not want to plot the graph without normalization
    # but the output of the binned array is needed for calculation
    # for location and confidence levels below
    pylab.close()

    fig = pylab.figure()
    # plot the normalized version (pdf) of the data
    pdf, bins, holder = pylab.hist(x, bins=N_bins, histtype='step',
                                   weights=prob, range=histrange,
                                   color='k', linewidth=2, normed=1)

    # Calculate the location and %confidence intervals
    # Since my location and confidence calculations can't take weighted
    # data I need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
        else:
            x_temp = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
            x_binned = numpy.concatenate((x_binned, x_temp))

    loc = biweightLoc(x_binned)
    ll_68, ul_68 = bcpcl(loc, x_binned, 1)
    ll_95, ul_95 = bcpcl(loc, x_binned, 2)

    # Create location and confidence interval line plots
    pylab.plot((loc, loc), (pylab.ylim()[0], pylab.ylim()[1]),
               '--k', linewidth=2, label='$C_{BI}$')
    pylab.plot((ll_68, ll_68), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#800000',
               label=r"68% $IC_{B_{BI}}$")
    pylab.plot((ul_68, ul_68), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#800000')

    pylab.plot((ll_95, ll_95), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#0000A0',
               label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95, ul_95), (pylab.ylim()[0], pylab.ylim()[1]), ':',
               linewidth=2, color='#0000A0')

    if x_label != None:
        pylab.xlabel(x_label, fontsize=15)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=15)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
        pylab.legend()
    if title != None:
        pylab.title(title, fontsize = 16)

    ### set font size
    # fontsize=14
    # ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)

    if save is True:
        filename = prefix + '_histplot1D'
        pylab.savefig(filename + '.pdf', bbox_inches='tight')

    if verbose is True:
        print '{0}, {1:0.4f}, {2:0.4f},'.format(prefix, loc, ll_68) + \
            '{0:0.4f}, {1:0.4f}, {2:0.4f}'.format(ul_68, ll_95, ul_95)

    return loc, ll_68, ul_68, ll_95, ul_95


# plots data after binning and weighting the data appropriately
def histplot1d(x, prefix, prob=None, norm=False, N_bins=100, histrange=None, x_lim=None, y_lim=None, x_label=None, y_label=None, legend=None):
    fig = pylab.figure()
    hist, binedges, tmp = pylab.hist(
        x, bins=N_bins, histtype='step', weights=prob, range=histrange, color='k', linewidth=2, normed=norm)

    # Calculate the location and %confidence intervals
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations

    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
        else:
            x_temp = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
            x_binned = numpy.concatenate((x_binned, x_temp))

    loc = biweightLoc(x_binned)
    ll_68, ul_68 = bcpcl(loc, x_binned, 1)
    ll_95, ul_95 = bcpcl(loc, x_binned, 2)
    # Create location and confidence interval line plots
    pylab.plot((loc, loc), (pylab.ylim()[0], pylab.ylim()[1]),
               '--k', linewidth=2, label='$C_{BI}$')
    pylab.plot((ll_68, ll_68), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#800000', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68, ul_68),
               (pylab.ylim()[0], pylab.ylim()[1]), '-.', linewidth=2, color='#800000')
    pylab.plot((ll_95, ll_95), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#0000A0', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95, ul_95),
               (pylab.ylim()[0], pylab.ylim()[1]), ':', linewidth=2, color='#0000A0')

    if x_label != None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
        pylab.legend()
    # fontsize=14
    #ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)

    filename = prefix+'_histplot1D'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(prefix, loc, ll_68, ul_68, ll_95, ul_95)

    return loc, ll_68, ul_68, ll_95, ul_95


# similar to histplot1d
def histplot1d_part(ax, x, prob=None, N_bins=100, histrange=None, x_lim=None, y_lim=None):
    '''
    This take the additional value of an array axes. for use with subplots
    '''
    hist, binedges, tmp = ax.hist(
        x, bins=N_bins, histtype='step', weights=prob, range=histrange, color='k', linewidth=1)

    # Calculate the location and %confidence intervals
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
        elif numpy.size(x_binned) == 0:
            x_binned = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
        else:
            x_temp = numpy.ones(hist[i])*(binedges[i]+binedges[i+1])/2
            x_binned = numpy.concatenate((x_binned, x_temp))
    loc = biweightLoc(x_binned)
    ll_68, ul_68 = bcpcl(loc, x_binned, 1)
    ll_95, ul_95 = bcpcl(loc, x_binned, 2)
    # Create location and confidence interval line plots
    ax.plot((loc, loc), (ax.get_ylim()[0], ax.get_ylim()[1]),
            '--k', linewidth=1, label='$C_{BI}$')
    ax.plot((ll_68, ll_68), (ax.get_ylim()[0], ax.get_ylim()[1]),
            '-.', linewidth=1, color='#800000', label='68% $IC_{B_{BI}}$')
    ax.plot((ul_68, ul_68),
            (ax.get_ylim()[0], ax.get_ylim()[1]), '-.', linewidth=1, color='#800000')
    ax.plot((ll_95, ll_95), (ax.get_ylim()[0], ax.get_ylim()[1]),
            ':', linewidth=1, color='#0000A0', label='95% $IC_{B_{BI}}$')
    ax.plot((ul_95, ul_95),
            (ax.get_ylim()[0], ax.get_ylim()[1]), ':', linewidth=1, color='#0000A0')

    if x_lim != None:
        ax.set_xlim(x_lim)
    if y_lim != None:
        ax.set_ylim(y_lim)
    return loc, ll_68, ul_68, ll_95, ul_95


# plot 2d histogram of 2 data arrays
def histplot2d(x, y, prefix, prob=None, N_bins=100, histrange=None, x_lim=None, y_lim=None, x_label=None, y_label=None, legend=None):
    '''
    Input:
    x = [1D array of N floats]
    y = [1D array of N floats]
    prefix = [string] prefix of output file
    prob = [None] or [1D array of N floats] weights to apply to each (x,y) pair
    N_bins = [integer] the number of bins in the x and y directions
    histrange = [None] or [array of floats: (x_min,x_max,y_min,y_max)] the range
        over which to perform the 2D histogram and estimate the confidence
        intervals
    x_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    y_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    x_label = [None] or [string] the plot's x-axis label
    y_label = [None] or [string] the plot's y-axis label
    legend = [None] or [True] whether to display a legend or not
    '''
    # Create the confidence interval plot
    if histrange == None:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]], weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins**2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum/h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum/h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1]+(xedges[1]-xedges[0])/2
    y = yedges[:-1]+(yedges[1]-yedges[0])/2
    X, Y = numpy.meshgrid(x, y)

    fig = pylab.figure()
    # Countours
    CS = pylab.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2))
    # imshow
    #im = pylab.imshow(H,cmap=pylab.cm.gray)
    pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r)

    if x_label != None:
        pylab.xlabel(x_label, fontsize=14)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=14)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
        # Dummy lines for legend
        # 800000 is for Maroon color - 68% 1 sigma
        # 0000A0 is for blue color - 95% confidence 3 sigma
        pylab.plot((0, 1), (0, 1), c='#800000', linewidth=2, label=('68%'))
        pylab.plot((0, 1), (0, 1), c='#0000A0', linewidth=2, label=('95%'))
        pylab.legend(scatterpoints=1)
    fontsize = 20
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    filename = prefix+'_histplot2d'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    return fig


# similar to histplot2d
def histplot2d_part(ax, x, y, prob=None, N_bins=100, histrange=None, x_lim=None, y_lim=None):
    '''
    This take the additional value of an array axes. for use with subplots
    Input:
    x = [1D array of N floats]
    y = [1D array of N floats]
    prefix = [string] prefix of output file
    prob = [None] or [1D array of N floats] weights to apply to each (x,y) pair
    N_bins = [integer] the number of bins in the x and y directions
    histrange = [None] or [array of floats: (x_min,x_max,y_min,y_max)] the range
        over which to perform the 2D histogram and estimate the confidence
        intervals
    x_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    y_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    x_label = [None] or [string] the plot's x-axis label
    y_label = [None] or [string] the plot's y-axis label
    legend = [None] or [True] whether to display a legend or not
    '''
    # Create the confidence interval plot
    if histrange == None:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]], weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins**2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum/h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum/h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1]+(xedges[1]-xedges[0])/2
    y = yedges[:-1]+(yedges[1]-yedges[0])/2
    X, Y = numpy.meshgrid(x, y)

    # Countours
    CS = ax.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2),
                    colors=((158/255., 202/255., 225/255.), (49/255., 130/255., 189/255.)))
    # imshow
    #im = ax.imshow(H,cmap=ax.cm.gray)
    ax.pcolor(X, Y, H, cmap=pylab.cm.gray)

    if x_lim != None:
        ax.set_xlim(x_lim)
    if y_lim != None:
        ax.set_ylim(y_lim)


# this is the one for generating the 2d plot in will 's paper?...
def histplot2dTSC(x, y, prefix, prob=None, N_bins=100, histrange=None, x_lim=None, y_lim=None, x_label=None, y_label=None, legend=None):
    '''
    Input:
    x = [1D array of N floats]
    y = [1D array of N floats]
    prefix = [string] prefix of output file
    prob = [None] or [1D array of N floats] weights to apply to each (x,y) pair
    N_bins = [integer] the number of bins in the x and y directions
    histrange = [None] or [array of floats: (x_min,x_max,y_min,y_max)] the range
        over which to perform the 2D histogram and estimate the confidence
        intervals
    x_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    y_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    x_label = [None] or [string] the plot's x-axis label
    y_label = [None] or [string] the plot's y-axis label
    legend = [None] or [True] whether to display a legend or not
    '''
    # Input calculated v and t parameters for other Dissociative Mergers
    v_bullet_analytic = 3400
    t_bullet_analytic = 0.218

    v_bullet_sf07 = 3400
    t_bullet_sf07 = 0.18

    v_macs = 2000
    t_macs = 0.255

    v_a520 = 2300
    t_a520 = 0.24

    v_pandora = 4045
    t_pandora = 0.162

    # Create the confidence interval plot
    if histrange == None:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]], weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins**2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum/h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum/h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1]+(xedges[1]-xedges[0])/2
    y = yedges[:-1]+(yedges[1]-yedges[0])/2
    X, Y = numpy.meshgrid(x, y)

    fig = pylab.figure()
    # Countours
    CS = pylab.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2))
    # imshow
    #im = pylab.imshow(H,cmap=pylab.cm.gray)
    pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r)

    # Data points for other dissociative mergers
    pylab.scatter(v_bullet_sf07, t_bullet_sf07, s=140,
                  c='k', marker='d', label="Bullet SF07")
    #pylab.scatter(v_macs,t_macs,s=140, c='0.4',markeredgecolor='0.4', marker='^',label='MACS J0025.4')
    # pylab.scatter(v_a520,t_a520,s=140,c='0.4',markeredgecolor='0.4',marker='o',label='A520')
    # pylab.scatter(v_pandora,t_pandora,s=140,c='0.4',markeredgecolor='0.4',marker='p',label='A2744')

    if x_label != None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
        # Dummy lines for legend
        pylab.plot((0, 1), (0, 1), c='#800000', linewidth=2, label=('68%'))
        pylab.plot((0, 1), (0, 1), c='#0000A0', linewidth=2, label=('95%'))
        pylab.legend(scatterpoints=1)
    fontsize = 20
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    filename = prefix+'_histplot2dTSC.pdf'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    return fig


# Plot the percentage differences between different % of data within
# Will's   2 million iterations
def plot_perdiff(perdiff, labels, title):
    fig = plt.figure()
    for i in range(5):
        plt.plot(range(1, len(perdiff[:, 1])+1),
                 perdiff[:, i], label=labels[i])
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc=4, bbox_to_anchor=(1.375, .5))
    minorLocator = MultipleLocator(2)

    plt.title(title)
    plt.grid()
    plt.grid(True, which='minor')
    plt.xlabel('% of iterations')
    plt.ylabel('Percent difference')
    ax.xaxis.set_minor_locator(minorLocator)
    plt.savefig(title+'.png', dpi=300, bbox_inches='tight')
    return fig


# Plot the histograms between the same set of data with
# and without applying prior on the same plot
def prior_diff(data1, data2, prefix, prob1=None, prob2=None, N_bins=100,
               histrange=None, x_lim=None, y_lim=None, x_label=None,
               y_label=None, legend=None):

    fig = pylab.figure()

    # bin data 1 first
    hist2, binedges2, tmp2 = pylab.hist(data2, bins=N_bins, histtype='step',
                                        weights=prob2, range=histrange, color='#ff0000', linewidth=2)
    hist1, binedges1, tmp1 = pylab.hist(data1, bins=N_bins, histtype='step',
                                        weights=prob1, range=histrange, color='#0000ff', linewidth=2)
    # print prefix+' hist1 array with len= ',len(hist1),' is '
    # print hist1
    # print prefix+' binedges1 is with len= ',len(binedges1),' is '
    # print binedges1
    # print prefix+' prob1 is with len= ',len(prob1),' is '
    # print prob1

    # bin data 2
    # print prefix+' hist2 array with len= ',len(hist2),' is '
    # print hist2
    # print prefix+' binedges2 is with len= ',len(binedges2),' is '
    # print binedges2
    # print prefix+' prob2 is with len= ',len(prob2),' is '
    # print prob2

    # Calculate the location and %confidence intervals for data 1
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned1 = numpy.ones(hist1[i])*(binedges1[i]+binedges1[i+1])/2
        else:
            x_temp1 = numpy.ones(hist1[i])*(binedges1[i]+binedges1[i+1])/2
            x_binned1 = numpy.concatenate((x_binned1, x_temp1))
    loc1 = biweightLoc(x_binned1)
    ll_68_1, ul_68_1 = bcpcl(loc1, x_binned1, 1)
    ll_95_1, ul_95_1 = bcpcl(loc1, x_binned1, 2)

    # Create location and confidence interval line plots
    pylab.plot((loc1, loc1), (pylab.ylim()[0], pylab.ylim()[1]),
               '--', linewidth=2, color='#6495ed', label='$C_{BI}$')
    pylab.plot((ll_68_1, ll_68_1), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#7fffd4', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68_1, ul_68_1),
               (pylab.ylim()[0], pylab.ylim()[1]), '-.', linewidth=2, color='#7fffd4')
    pylab.plot((ll_95_1, ll_95_1), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#87ceeb', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95_1, ul_95_1),
               (pylab.ylim()[0], pylab.ylim()[1]), ':', linewidth=2, color='#87ceeb')

    # Calculate the location and %confidence intervals for data 2
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned2 = numpy.ones(hist2[i])*(binedges2[i]+binedges2[i+1])/2
        else:
            x_temp2 = numpy.ones(hist2[i])*(binedges2[i]+binedges2[i+1])/2
            x_binned2 = numpy.concatenate((x_binned2, x_temp2))
    loc2 = biweightLoc(x_binned2)
    ll_68_2, ul_68_2 = bcpcl(loc2, x_binned2, 1)
    ll_95_2, ul_95_2 = bcpcl(loc2, x_binned2, 2)

    # Create location and confidence interval line plots
    pylab.plot((loc2, loc2), (pylab.ylim()[0], pylab.ylim()[1]),
               '--', linewidth=2, color='#ff4500', label='$C_{BI}$')
    pylab.plot((ll_68_2, ll_68_2), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#ff8c00', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68_2, ul_68_2),
               (pylab.ylim()[0], pylab.ylim()[1]), '-.', linewidth=2, color='#ff8c00')
    pylab.plot((ll_95_2, ll_95_2), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#ffa500', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95_2, ul_95_2),
               (pylab.ylim()[0], pylab.ylim()[1]), ':', linewidth=2, color='#ffa500')

    # create labels for the plots
    if x_label != None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
        pylab.legend()
    # fontsize=14
    #ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)

    filename = prefix+'_prior_diff'
    pylab.savefig(title+'.png', dpi=300, bbox_inches='tight')

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(prefix, loc1, ll_68_1, ul_68_1, ll_95_1, ul_95_1)

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(prefix, loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2)

    return loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2

# Plot the pdf between the same set of data with
# and without applying prior on the same plot


def prior_diff_pdf(data1, data2, prefix, prob1=None, prob2=None, N_bins=100,
                   histrange=None, x_lim=None, y_lim=None, x_label=None,
                   y_label=None, legend=None):

    # bin data 2 first
    hist2, binedges2, tmp2 = pylab.hist(data2, bins=N_bins, histtype='step',
                                        weights=prob2, range=histrange, color='#ff0000', linewidth=2)
    # bin data 1
    hist1, binedges1, tmp1 = pylab.hist(data1, bins=N_bins, histtype='step',
                                        weights=prob1, range=histrange, color='#0000ff', linewidth=2)
    pylab.close()

    fig = pylab.figure()
    # plot the pdf for the data
    pdf2, histbin2, tmp_pdf2 = pylab.hist(data2, bins=N_bins, normed=1,
                                          histtype='step', weights=prob2, range=histrange, color='#ff0000',
                                          linewidth=2)

    pdf1, histbin1, tmp_pdf1 = pylab.hist(data1, bins=N_bins, normed=1,
                                          histtype='step', weights=prob1, range=histrange, color='#0000ff',
                                          linewidth=2)

    # Calculate the location and %confidence intervals for data 1
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned1 = numpy.ones(histbin1[i]) *\
                (binedges1[i]+binedges1[i+1])/2
        else:
            x_temp1 = numpy.ones(hist1[i])*(binedges1[i]+binedges1[i+1])/2
            x_binned1 = numpy.concatenate((x_binned1, x_temp1))
    loc1 = biweightLoc(x_binned1)
    ll_68_1, ul_68_1 = bcpcl(loc1, x_binned1, 1)
    ll_95_1, ul_95_1 = bcpcl(loc1, x_binned1, 2)

    # adjust the max ylim so it does not look weird
    ylim_max = pdf2.max()*1.2

    # Create location and confidence interval line plots
    pylab.plot((loc1, loc1), (pylab.ylim()[0], ylim_max), '--', linewidth=2,
               color='#6495ed', label='$C_{BI}$')
    pylab.plot((ll_68_1, ll_68_1), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#7fffd4', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68_1, ul_68_1), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#7fffd4')
    pylab.plot((ll_95_1, ll_95_1), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#87ceeb', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95_1, ul_95_1), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#87ceeb')

    # Calculate the location and %confidence intervals for data 2
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned2 = numpy.ones(hist2[i])*(binedges2[i]+binedges2[i+1])/2
        else:
            x_temp2 = numpy.ones(hist2[i])*(binedges2[i]+binedges2[i+1])/2
            x_binned2 = numpy.concatenate((x_binned2, x_temp2))
    loc2 = biweightLoc(x_binned2)
    ll_68_2, ul_68_2 = bcpcl(loc2, x_binned2, 1)
    ll_95_2, ul_95_2 = bcpcl(loc2, x_binned2, 2)

    # Create location and confidence interval line plots
    # if y_lim == None:
    # else:
    #    ylim_max = y_lim[1]
    #    ylim_min = y_lim[0]
    # if x_lim == None:
    #    xlim_max = pylab.xlim()[1]
    #    xlim_min = pylab.xlim()[0]
    # else:
    #    xlim_max = x_lim[1]
    #    xlim_min = x_lim[0]
    pylab.plot((loc2, loc2), (pylab.ylim()[0], ylim_max), '--', linewidth=2,
               color='#ff4500', label='$C_{BI}$')
    pylab.plot((ll_68_2, ll_68_2), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#ff8c00', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68_2, ul_68_2), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#ff8c00')
    pylab.plot((ll_95_2, ll_95_2), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#ffa500', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95_2, ul_95_2), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#ffa500')

    # create labels for the plots
    if x_label != None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
        pylab.legend()
    # fontsize=14
    #ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)
    pylab.ylim(0, pdf2.max()*1.2)

    filename = prefix+'_prior_diff'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(prefix, loc1, ll_68_1, ul_68_1, ll_95_1, ul_95_1)

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(prefix, loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2)

    return loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2

# My function for plotting 2 sets of 2d contour


def histplot2d_2contour(
    x1, y1, x, y, prefix, prob1=None, prob=None, N_bins=100,
    histrange=None, x_lim=None, y_lim=None, x_label=None, y_label=None,
        legend=None):
     # Create the confidence interval plot
    if histrange == None:
        if prob1 != None:
            H1, xedges1, yedges1 = numpy.histogram2d(x1, y1, bins=N_bins,
                                                     weights=prob1)
        elif prob1 == None:
            H1, xedges1, yedges1 = numpy.histogram2d(x1, y1, bins=N_bins)
    else:
        if prob1 != None:
            H1, xedges1, yedges1 = numpy.histogram2d(x1, y1, bins=N_bins,
                                                     range=[
                                                         [histrange[0],
                                                             histrange[1]],
                                                         [histrange[2], histrange[3]]], weights=prob1)
        elif prob == None:
            H1, xedges1, yedges1 = numpy.histogram2d(x1, y1, bins=N_bins,
                                                     range=[
                                                         [histrange[0],
                                                             histrange[1]],
                                                         [histrange[2], histrange[3]]])
    H1 = numpy.transpose(H1)
    # Flatten H
    h = numpy.reshape(H1, (N_bins**2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum/h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum/h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]

    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x1 = xedges1[:-1]+(xedges1[1]-xedges1[0])/2
    y1 = yedges1[:-1]+(yedges1[1]-yedges1[0])/2
    X1, Y1 = numpy.meshgrid(x1, y1)

    fig = pylab.figure()
    # Countours
    CS = pylab.contour(X1, Y1, H1, (h_2sigma, h_1sigma),
                       linewidths=(2, 2), colors=('#a4a4a4', '#6e6e6e'))
    # imshow
    #im = pylab.imshow(H1,cmap=pylab.cm.gray)
    # pylab.pcolor(X1,Y1,H1,cmap=pylab.cm.white)

    if x_label != None:
        pylab.xlabel(x_label, fontsize=14)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=14)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    # if legend != None:
    # Dummy lines for legend
        # 800000 is for light gray - 68% 1 sigma
        # 0000A0 is for whitesmoke - 95% confidence 3 sigma
    # pylab.plot((0,1),(0,1),c='#f5fffa',linewidth=2,label=('68%'))
    # pylab.plot((0,1),(0,1),c='#d3d3d3',linewidth=2,label=('95%'))
        # pylab.legend(scatterpoints=1)
    fontsize = 20
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    # SECOND contour
    # Create the confidence interval plot for the second sets of contour
    if histrange == None:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob != None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]], weights=prob)
        elif prob == None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]], [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins**2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum/h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum/h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1]+(xedges[1]-xedges[0])/2
    y = yedges[:-1]+(yedges[1]-yedges[0])/2
    X, Y = numpy.meshgrid(x, y)

    # Coutours
    CS = pylab.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2))
    # imshow
    #im = pylab.imshow(H,cmap=pylab.cm.Blues)
    pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r)

    if x_label != None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label != None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim != None:
        pylab.xlim(x_lim)
    if y_lim != None:
        pylab.ylim(y_lim)
    if legend != None:
    # Dummy lines for legend
        # 800000 is for Maroon color - 68% 1 sigma
        # 0000A0 is for blue color - 95% confidence 3 sigma
        pylab.plot((0, 1), (0, 1), c='#87cefa', linewidth=2, label=('68%'))
        pylab.plot((0, 1), (0, 1), c='#6495ed', linewidth=2, label=('95%'))
        pylab.legend(scatterpoints=1)
    fontsize = 20
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    filename = prefix+'_histplot2d_combine2'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    return fig


# Plot the percentage differences between different % of data within
# Will's   2 million iterations
def plot_perdiff(perdiff, labels, title):
    fig = plt.figure()
    for i in range(5):
        plt.plot(range(1, len(perdiff[:, 1])+1),
                 perdiff[:, i], label=labels[i])
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc=4, bbox_to_anchor=(1.375, .5))
    minorLocator = MultipleLocator(2)

    plt.title(title)
    plt.grid()
    plt.grid(True, which='minor')
    plt.xlabel('% of iterations')
    plt.ylabel('Percent difference')
    ax.xaxis.set_minor_locator(minorLocator)
    plt.savefig(title+'.png', dpi=300, bbox_inches='tight')
    return fig
