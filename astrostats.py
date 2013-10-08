'''
astrostats.py is written Will Dawson and contains various statistical functions
that I have found useful (e.g. the Beers et al. 1990 paper).
-- On 5/13/2012 I corrected an error in the MAD calculations, it was MAD = numpy.median(numpy.abs(z))
'''
from __future__ import division
import numpy
from scipy.stats import norm
from scipy.special import erf
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

"""
Copyright (c) 2012, William A. Dawson
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the UC Davis nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL WILLIAM A. DAWSON BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
