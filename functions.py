import numpy
from scipy.interpolate import interp1d
from scipy.stats import norm
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

def wave(x, amplitude, frequency, phi, costant): 
    return amplitude * numpy.cos(frequency * x + phi) + costant

def increasing_wave(x, amplitude, frequency, phi, costant, m): 
    return amplitude * numpy.cos(frequency * x + phi) + costant + m * x

def gauss(x, norm, mean, sigma): 
    return (norm/(sigma * numpy.sqrt(2 * numpy.pi))) * numpy.exp(-0.5 * ((x - mean)/sigma )**2)
  
def two_gauss(x, a, norm, mean1, sigma1, mean2, sigma2):
    return a * gauss(x, norm, mean1, sigma1) + (1.-a) * gauss(x, norm, mean2, sigma2)

def exponential(x, a, m, costant): 
    return  (a/m ) * numpy.exp(-x / m) + costant

def two_expo(x, norm, fraction, m_short, m_long, costant): 
    return  norm * (fraction * numpy.exp(- x / m_short) + (1. - fraction) * numpy.exp(- x / m_long) ) + costant      
#def two_expo(x, a1, m_short, a2, m_long, costant):
 #   return a1 * numpy.exp(-x / m_short) + a2 * numpy.exp(-x/m_long) + costant
    
    
def two_expo_gauss(x, norm, fraction, m_short, m_long, costant, gauss_norm, gauss_mean, gauss_sigma):     
    return two_expo(x, norm, fraction, m_short, m_long, costant) + gauss(x, gauss_norm, gauss_mean, gauss_sigma) 

def two_expo_integral(bin_center, norm, fraction, m_short, m_long, costant):
    t_sup = bin_center + 0.5 * (bin_center[1] - bin_center[0])
    t_inf = bin_center - 0.5 * (bin_center[1] - bin_center[0])
    return norm * 1000 * (fraction * m_short * (numpy.exp(-t_inf/m_short) - numpy.exp(-t_sup/m_short)) + (1.- fraction) * m_long * (numpy.exp(-t_inf/m_long) - numpy.exp(-t_sup/m_long))) + costant * (t_sup - t_inf)


def gauss_log_likelihood(x, y, sigma, model, *model_params):
    y_pred = model(x, *model_params)
    #print('ypred', y_pred)
    #print('norm.pdf(y, loc=y_pred, scale=sigma)', norm.pdf(y, loc=y_pred, scale=sigma))
    return numpy.sum(numpy.log(norm.pdf(y, loc=y_pred, scale=sigma)))

def ll_ratio_test_stat(loglh_alt, loglh_null):
    return -2 * (loglh_null - loglh_alt)    
    

def poisson_log_likelihood(bins_center, n, model): 
    def likelihood(model_params):
        y_pred = model(bins_center, *tuple(model_params))
        return numpy.sum(y_pred - n + n * (numpy.log(n) - numpy.log(y_pred)))
    return likelihood

def two_expo_integral_jacobian(bin_center, n):
    bin_width = bin_center[1] - bin_center[0]
    print('bin_width = ', bin_width)
    def jacobian(par):
        norm, fraction, m_short, m_long, costant = par[0], par[1], par[2], par[3], par[4]
        pred = two_expo_integral(bin_center, norm, fraction, m_short, m_long, costant)
        t_sup = bin_center + 0.5 * bin_width
        t_inf = bin_center - 0.5 * bin_width
        exp_sup_short = numpy.exp(-t_sup/m_short)
        exp_inf_short = numpy.exp(-t_inf/m_short)
        exp_sup_long =  numpy.exp(-t_sup/m_long)
        exp_inf_long =  numpy.exp(-t_inf/m_long)
        
        first = 1000 * (fraction * m_short * (exp_inf_short - exp_sup_short) + (1.- fraction) * m_long * (exp_inf_long - exp_sup_long))
        second = 1000 * (norm * m_short * (exp_inf_short - exp_sup_short) - norm * m_long * (exp_inf_long - exp_sup_long))
        third = 1000 * (norm * fraction * (exp_inf_short - exp_sup_short) + norm * fraction * (t_inf/m_short * exp_inf_short - t_sup/m_short * exp_sup_short ))
        fourth = 1000 * (norm * (1. - fraction) * (exp_inf_long- exp_sup_long) + norm * (1. - fraction) * (t_inf/m_long * exp_inf_long - t_sup/m_long * exp_sup_long))
        fifth = (t_sup - t_inf)
        return numpy.sum((1. - n/pred) * numpy.array([first, second, third, fourth, fifth]), axis=1)
    return jacobian






























