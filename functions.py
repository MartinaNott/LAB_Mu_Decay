import numpy
from scipy.interpolate import interp1d
import scipy.special

import constants

#Distributions for the generation of the events:
def dist_theta(theta):
  return (numpy.cos(theta))**2
  
def power_law(x, index, x0=1.):
    return (x/x0)**index

def positive_erf(x, mean, sigma):
    return 0.5 + 0.5 * scipy.special.erf((x - mean)/sigma)

def spectrum(x, index, x_drop, drop_width = 1.):
    return positive_erf(x, x_drop, drop_width)  * power_law(x, index)


def electron_spectrum(E_electron ): 
    x = 2. * E_electron/ constants.MUON_MASS
    return (3. - 2. * x ) * x**2


def muon_spectrum(E_muon): 
    return 0.14 * E_muon**2.7 * 0.001 * (1/( 1 + 1.1 * E_muon/0.115 ) + 0.054/(1 + 1.1 * E_muon/0.850))
    
    
#Fit functions:   
def line(x, m , q):
    return m * x +q

def proportional(x, m ):
    return m * x

def costant(x,  q): 
    q = numpy.ones(len(x))*q
    return q

def gauss(x, norm, mean, sigma): 
    return (norm/(sigma * numpy.sqrt(2 * numpy.pi))) * numpy.exp(-0.5 * ((x - mean)/sigma )**2)
  
def two_gauss(x, a, norm, mean1, sigma1, mean2, sigma2):
    return a * gauss(x, norm, mean1, sigma1) + (1.-a) * gauss(x, norm, mean2, sigma2)

def exponential(x, a, m, costant): 
    return  (a/m ) * numpy.exp(-x / m) + costant

def two_expo(x, norm, fraction, m_short, m_long, costant): 
    return  norm * (fraction * numpy.exp(- x / m_short) + (1. - fraction) * numpy.exp(- x / m_long) ) + costant      
    
    
def two_expo_gauss(x, norm, fraction, m_short, m_long, costant, gauss_norm, gauss_mean, gauss_sigma):     
    return two_expo(x, norm, fraction, m_short, m_long, costant) + gauss(x, gauss_norm, gauss_mean, gauss_sigma) 


def do_expo_integral():
    def expo_integral(bin_center, norm, fraction, m_short, m_long, costant):
        t_sup = bin_center + 0.5 * (bin_center[1] - bin_center[0])
        t_inf = bin_center - 0.5 * (bin_center[1] - bin_center[0])
        return - norm * (fraction/m_short * (numpy.exp(-t_sup/m_short) - numpy.exp(-t_inf/m_short)) + (1.- fraction)/m_long * (numpy.exp(-t_sup/m_long) - numpy.exp(-t_inf/m_long))) + costant * (t_sup - t_inf)
    return expo_integral

