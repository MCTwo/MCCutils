import numpy
import scipy.integrate
#constants
G = 6.673*10**-11 # m**3/kg/s**2
c = 3*10**5 #Units km/s
#conversions
secinGyr = 31556926*10**9 #seconds in a gigayear
kminMpc = 3.08568025*10**19 # km in a Megaparsec
kginMsun = 1.98892*10**30



def chi(x,Om=.3,Ol=.7):
    # proportional to root(1/H(t))
    return 1/numpy.sqrt(Ol+Om*(1+x)**3+(1-Ol-Om)*(1+x)**2)

def Dcom(z,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: Dcom(z,h=0.7,Om=0.3,Ol=0.7)

    Returns the radial comoving distance (in units of Mpc)
    """
    f0 = scipy.integrate.quad(lambda x: chi(x,Om,Ol),0,z)[0]
    da = c/(100.0*h)*f0
    return da

def Da(z,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: Da(z,h=0.7,Om=0.3,Ol=0.7)

    Returns the cosmological angular diameter distance (in units of Mpc)
    """
    f0 = scipy.integrate.quad(lambda x: chi(x,Om,Ol),0,z)[0]
    da = c/(100.0*h)/(1+z)*f0
    return da

def Dl(z,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: Dl(z,h=0.7,Om=0.3,Ol=0.7)

    Returns the cosmological luminosity distance (in units of Mpc)
    """
    return Da(z,h,Om,Ol)*(1+z)**2
    
def ProjectedLength(z,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: ProjectedLength(z,h=0.7,Om=0.3,Ol=0.7)

    Returns the projected length (in units of Mpc/arcmin)
    Assumes small angles (i.e.: Da*Theta = Projected Length, i.e. sin(Theta)=Theta)
    """
    return Da(z,h,Om,Ol)/3437.75 #given that 1 radian = 3437.75 arcminutes

def lensgeo(zl,zs,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: lensgeo(zl,zs,h=0.7,Om=0.3,Ol=0.7)

    Returns a dictionary containing the cosmological angular diameter distances
    (in units of Mpc) used in gravitational lensing analyis. As well as the
    critical surface mass density (in units kg/m**2).

    {'Dl':###,'Ds':###,'Dls':###,'sigcr':###}
    Dl = angular diameter distance to lens
    Ds = angular diameter distance to source
    Dls = angular diameter distance from lens to source
    sigcr = Critical surface mass density
    """
    if zl >= zs:
        print 'Error: Lens redshift must be smaller than source redshift.'
    else:
        f0l = scipy.integrate.quad(lambda x: chi(x,Om,Ol),0,zl)[0]
        f0s = scipy.integrate.quad(lambda x: chi(x,Om,Ol),0,zs)[0]
        
        ds = c/(100.0*h)/(1+zs)*f0s # distance to source
        dl = c/(100.0*h)/(1+zl)*f0l # distance to lens
        dls = c/(100.0*h)/(1+zs)*(f0s-f0l) # distance between lens and source
        
        sigcr = c**2/(4*numpy.pi*G)*ds/(dl*dls)*1000/kminMpc 

        return {'Dl':dl,'Ds':ds,'Dls':dls,'sigcr':sigcr}

def H(z,h=0.7,Om=0.3,Ol=0.7,Or=0,Ok=0):
    """
    Usage: H(z,h=0.7,Om=0.3,Ol=0.7,Or=0,Ok=0)

    Returns the hubble constant for the reshift given in Units (km/s)/Mpc
    """
    return numpy.sqrt((100*h)**2*(Ol+Om*(1+z)**3+Ok*(1+z)**2+Or*(1+z)**4))

def age(z,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: age(a,h=0.7,Om=0.3,Ol=0.7)

    Returns the approximate age of the universe at the given redshift
    Units Gyr
    """
    a = 1/(1.0+z)
    def Ht(x,Om=.3,Ol=.7):
        # proportional to root(1/H(t))
        return 1/numpy.sqrt(Ol*x**2+Om/x+(1-Ol-Om))
    h_Gyr = (100*h)/(3.08568*10**19)*3.1557*10**7*10**9
    t = 1/h_Gyr*scipy.integrate.quad(lambda x: Ht(x,Om,Ol),0,a)[0]
    return t

def lookbacktime(z,h=0.7,Om=0.3,Ol=0.7):
    """
    Usage: age(z,h=0.7,Om=0.3,Ol=0.7)

    Returns the approximate look-back time to the given redshift
    Units Gyr
    """
    return age(0,h,Om,Ol) - age(z,h,Om,Ol)
    
## These are not producing the correct answer.
##def age(z,h=0.7,Om=0.3,Ol=0.7,Or=0):
##    """
##    Usage:
##        age(z,h=0.7,Om=0.3,Ol=0.7,Or=0)
##
##    Returns the age of the universe in Gyr
##    """
##    return 1/(100*h)*kminMpc/secinGyr*scipy.integrate.quad(lambda x: (Or*(1+x)**6 + Om*(1+x)**5+(1-Om-Ol)*(1+x)**4+Ol*(1+z)**2)**(-1/2),z,scipy.integrate.Inf)[0]
##
##def lookback(z,h=0.7,Om=0.3,Ol=0.7,Or=0):
##    """
##    Usage:
##        age(z,h=0.7,Om=0.3,Ol=0.7,Or=0)
##
##    Returns the age of the universe in Gyr
##    """
##    return 1/(100*h)*kminMpc/secinGyr*scipy.integrate.quad(lambda x: (Or*(1+x)**6 + Om*(1+x)**5+(1-Om-Ol)*(1+x)**4+Ol*(1+z)**2)**(-1/2),0,z)[0]


def rhoCrit(z,h=0.7,Om=0.3,Ol=0.7,Or=0):
    """
    Usage:
        rhoCrit(z,h=0.7,Om=0.3,Ol=0.7,Or=0)

    Returns the critical density for of the universe for redshift z (in units kg/m**3)
    """
    return 3*(H(z,h,Om,Ol,Or)/kminMpc)**2/(8*numpy.pi*G)

def r200(m200,z,h=0.7,Om=0.3,Ol=0.7,Or=0):
    '''
    Usage:
        r200(m200,z,h=0.7,Om=0.3,Ol=0.7,Or=0)
    m200 is approximatly the viral mass of the object in solar mass units.
    z is the redshift of the mass.

    Returns the R200 radius for the mass in Mpc.
    '''
    rho = rhoCrit(z,h,Om,Ol,Or)/kginMsun*(1000*kminMpc)**3
    return (3*m200/(4.*numpy.pi*200.*rho))**(1/3.)

def v200(m200,z,h=0.7,Om=0.3,Ol=0.7,Or=0,Ok=0):
    '''
    Usage:
        v200(m200,z,h=0.7,Om=0.3,Ol=0.7,Or=0,Ok=0)
    m200 is the viral mass of the object in solar mass units.
    z is the redshift of the mass.

    Returns the virial velocity of the mass in km/s.
    '''
    return (m200*10*G*H(z,h,Om,Ol,Or,Ok)/kminMpc*kginMsun)**(1/3.)/1000.

def vdisp(m200,z,h=0.7,Om=0.3,Ol=0.7,Or=0,Ok=0):
    '''
    Usage:
        vdisp(m200,z,Rg=1,h=0.7,Om=0.3,Ol=0.7,Or=0,Ok=0)
    m200 is the M_{200} mass of the object in solar mass units.
    z is the redshift of the mass.

    Uses the relation in: Evrard, A.E. et al., 2008. Virial Scaling of Massive Dark Matter Halos: Why Clusters Prefer a High Normalization Cosmology. The Astrophysical Journal, 672(1), pp.122

    Returns the velocity dispersion of the mass in km/s.
    '''
    return 10**(numpy.log10(1080)+0.352*numpy.log10(H(z,h,Om,Ol,Or,Ok)/100.*m200/10^15))

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
    * Neither the name of the University of California, Davis nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
