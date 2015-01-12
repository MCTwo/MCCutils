'''
astrostats.py is written Will Dawson and contains various statistical functions
that I have found useful (e.g. the Beers et al. 1990 paper).
-- On 5/13/2012 I corrected an error in the MAD calculations, it was MAD = numpy.median(numpy.abs(z))
'''
from __future__ import division
import numpy
from scipy.stats import norm
from scipy.special import erf
import sys

def biweightLoc(z,c=6):
    '''
    Biweight statistic Location (similar to the mean).
    '''
    M = numpy.median(z)
    MAD = numpy.median(numpy.abs(z-M))
    u = (z-M)/(c*MAD)
    mask_u = numpy.abs(u) < 1
    z = z[mask_u]
    u = u[mask_u]
    Cbi = M + numpy.inner(z-M,(1-u**2)**2)/numpy.sum((1-u**2)**2)
    return Cbi

def biweightScale(z,c=9):
    '''
    Biweight statistic Scale (similar to the standard deviation).
    '''
    n = numpy.size(z)
    M = numpy.median(z)
    MAD = numpy.median(numpy.abs(z-M))
    u = (z-M)/(c*MAD)
    mask_u = numpy.abs(u) < 1
    z = z[mask_u]
    u = u[mask_u]
    Sbi = n**(0.5)*numpy.inner((z-M)**2,(1-u**2)**4)**(0.5)/numpy.abs(numpy.inner(1-u**2,1-5*u**2))
    return Sbi

def bcpcl(T,T_p,N_sigma):
    '''
    Calculates the bias corrected percent confidence limits.
    -- Suppose that we have observed data (y1, y2, ..., yn) and use it to estimate a population parameter Q (e.g. Q could be the true mean of the entire population).
    -- T is a statistic that estimates Q. For example T could be an estimate of the true mean by calculating the mean of  (y1, y2, ..., yn).
    -- Suppose that we create m bootstrap samples (y_p_1j, y_p_2j, ...,j_p_nj) from observed sample  (y1, y2, ..., yn), where j is the jth bootstrap sample.
    -- Then T_p_j is the jth bootstrap observation of T.  For example this could be the mean of (y_p_1j, y_p_2j, ...,j_p_nj).
    
    T = [float] e.g. biweight Location for (y1, y2, ..., yn)
    T_p = [vector array] biwieght Locations for the bootstrap samples
    N_sigma = the number of sigma to report the confidence limits for
        e.g. for 95% confidence limits N_sigma=2
    Return (lower, upper) confidence limits
    '''
    #Percentile confidence interval is defined as 100%(1-a), thus for 1sigma a=0.32
    a = 1-erf(N_sigma/numpy.sqrt(2))
    #order the bootstrap sample values smallest to largest
    index = numpy.argsort(T_p)
    T_p = T_p[index]
    #Number of bootstrap samples
    m = numpy.size(T_p)        
    #Calculate the bias correction term
    mask = T_p < T
    z_0 = norm.ppf(numpy.sum(mask)/m)
    #Calculate the a1 and a2 values
    a1 = norm.cdf(2*z_0+norm.ppf(a/2))
    a2 = norm.cdf(2*z_0+norm.ppf(1-a/2))
    #Calculate the lower and upper indicies of lower and upper confidence intervals
    id_L = numpy.int(m*a1)-1
    id_U = numpy.int(m*a2)
    #Find the lower an upper confidence values
    T_L = T_p[id_L]
    T_U = T_p[id_U]
    return T_L, T_U

def weightedrand(weights,size=1):
    '''
    Given a list of weights [w_0, w_1, ..., w_n-1] corresponding to values in an
    array [a_0, a_1, ..., a_n-1] it will return an index of a random draw
    with probability proportional to weights [w_0, w_1, ..., w_n-1].
    
    Input:
    weights = 1D array of length n
    
    Output:
    index = array index of randomly drawn value with probability 
      proportional to weights
    '''
    # check to make sure that there are no negative weights
    if numpy.sum(weights<0) > 0:
        print 'weightedrandom: error, negative weights are not allowed, exiting'
        sys.exit()
    # determine the normalizes cumulative sum of the weight array    
    w_norm = weights/numpy.sum(weights)
    w_cumsum = numpy.cumsum(w_norm)
    # build the sampled array
    if size == 1:
        # randomly draw a number from 0-1
        rand = numpy.random.rand()
        # find where this random number intersects the cummulative weight function
        diff = w_cumsum-rand
        index = (1/diff).argmin()
    else:
        index = numpy.zeros(size,dtype=numpy.uint64)
        # randomly draw a number from 0-1
        rand = numpy.random.rand(size)
        for i in range(size):
            diff = w_cumsum-rand[i]
            index[i] = (1/diff).argmin()
    return index
