import sys
sys.path.insert(1, '/home/elo/ele/LAB_Mu_Decay')

import matplotlib.pyplot as plt
import numpy
import argparse 
import datetime
import time
import scipy.integrate as integrate
from itertools import product
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import minimize, check_grad, approx_fprime

import event_generator_functions
import geometry
import functions
import constants
import plot_functions
import utilities

def mu_decay(num_events, pdf, t_min, t_max, *args): 
  dt = numpy.linspace(t_min, t_max, 1000)  
  cdf_y  = numpy.full(len(dt), 0.)
  for i in range(len(dt)):
    y, rest = quad(pdf, dt[0], dt[i], args = args)
    cdf_y[i] = y       
  cdf_y, unique_indices = numpy.unique(cdf_y, return_index=True)
  cdf_y = cdf_y/cdf_y[-1]
  dt = dt[unique_indices]  
  ppf_spline = interp1d(cdf_y, dt)         
  x = numpy.random.uniform(0., 1., num_events)
  delta_t_decay = ppf_spline(x)  
  return delta_t_decay



def likelihood_fit(model, time_diff=None, bin_center=None, n=None, jacobian = None,  x0=None, bounds=None, title='', legend='', fit_min = -numpy.inf, range_hist = (0., 20), n_bins = 100):

    plt.subplot(2,1,1)
    if time_diff is not None: 
        bin_center, n, dn = plot_functions.plot_histogram(time_diff, "time [$\mu$s]", "", n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True)  
    else: 
        dn = numpy.sqrt(n)
        plot_functions.scatter_plot(bin_center, n, "time [$\mu$s]", "", dy = dn , title = title, fmt = '.b')    
    
    mask = bin_center >  fit_min   
    bin_center = bin_center[mask]
    n = n[mask]
    minus_two_ll = functions.poisson_log_likelihood(bin_center, n, model)
    if jacobian is not None:
        jac = jacobian(bin_center, n)
    else:
        jac = None
    res = minimize(minus_two_ll, x0 = x0, bounds = bounds , method = 'BFGS', jac = jac, options={'disp' : True})
    opt = res.x 
    pcov = res.hess_inv
    print('OPT_likelihood 2expo:', opt)
    print(numpy.sqrt(numpy.diagonal(pcov)/len(n)))
    plot_functions.scatter_plot(bin_center, model(bin_center, *opt), "time [$\mu$s]", "", fmt='-')   
    plt.subplot(2,1,2)
    residuals = n - model(bin_center, *opt)
    plot_functions.scatter_plot(bin_center, residuals, "time [$\mu$s]", "res", fmt='.')    
    return minus_two_ll, bin_center, n, dn

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('--number_events', '-n', default=1000, type=int, help='')
options_parser.add_argument('--output_file', '-o', default='None', type=str, help='')


if __name__ == '__main__' :   
    numpy.random.seed(3)
        
    start_time = time.time()   
    options = vars(options_parser.parse_args())  
    N = options['number_events']
    output_file_events = options['output_file']

    param_names_2exp = ['norm', 'fraction', 'm_short', 'm_long', 'costant']
    param_units_2exp = ['$\mu ^-1$s', '', '$\mu$s', '$\mu$s', '$\mu ^-1$s']
    
    #norm, fraction, m_short, m_long, costant
    p0_gen = [1.586, 0.50054369 , 0.088712, 1.97491405, 0.007]
    n_bins = 150
    fit_min = 0.000

    """
    data_file = 'data/run11_mudecay.dat', 'data/run9_mudecay.dat' #10512 eventi up e 1657 eventi down 
    data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file])
    ch = data[0, :]
    time = data[1, :]
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch, time, 5, 7)   
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch, time, 5, 8)   
    time_diff = numpy.concatenate((time_diff_up, time_diff_down))
    
    mask_start = ch == 5
    a = numpy.ones(len(ch))
    a = a[mask_start]
    print('start down', len(a)) 
    bins = numpy.linspace(0., 15., n_bins + 1)
    x = 0.5 * (bins[1:] + bins[:-1])
    """
    norm_noise = 2.39 * 3291430/85227 #3291430/85227
    tau_noise = 0.106
    cost_noise = 1.2
    N_noise = int(100 * 3291430/85227)
    N_events = N - N_noise
    p0_noise_gen = [norm_noise, tau_noise, cost_noise]
    dt = mu_decay(N_events, functions.two_expo, 0., 100., *p0_gen)
    dt_noise = mu_decay(N_noise, functions.exponential, 0., 100., *p0_noise_gen)
    dt = numpy.concatenate((dt, dt_noise))
    #dt_noise = mu_decay(230, functions.two_expo, 0.075, 0.2250, *p0_gen)
    #dt = numpy.concatenate((dt, dt_noise))

        
    plt.figure()
    p0 = [6.2, 0.50054369 , 0.0882712, 1.97197491405, 0.7]
    minus_two_ll, bin_center, n, dn = likelihood_fit( model = functions.two_expo_integral,
                    time_diff= dt, jacobian = functions.two_expo_integral_jacobian, 
                    n_bins = n_bins, x0 = p0, bounds = None, fit_min = fit_min,  
                    range_hist = (0., 15.), title ='MC, con noise')

    plt.subplot(2,1,1)
    p0 = [7, 0.50054369 , 0.0882712, 2.197491405, 0.7]
    bounds =  (0.0, 0.01, 0.02, 1.5, 0.), (numpy.inf, 1., 0.300, 5., 1000)
    opt, pcov = plot_functions.do_fit(bin_center, n, dn, param_names=param_names_2exp, param_units=param_units_2exp, fit_function = functions.two_expo_integral, p0 = p0, show=True, bounds = bounds, draw_on_points=True)
    print(pcov)

    """
    plt.figure()
    n = n - functions.exponential(bin_center, norm_noise, tau_noise, cost_noise)
    minus_two_ll, bin_center, n, dn = likelihood_fit( model = functions.two_expo_integral,
                    bin_center=bin_center, n = n, jacobian = 
                    functions.two_expo_integral_jacobian, 
                    n_bins = n_bins, x0 = p0, bounds = None, fit_min = 0.10,  
                    range_hist = (0., 20.), title ='MC, senza noise')
    
    plt.figure()
    plt.title('confronto')
    plot_functions.plot_histogram(dt, '', '', n_bins = n_bins, range = (0., 15.), title = '', legend = 'mc', fmt = '.g', as_scatter = True)    
    plot_functions.scatter_plot(x, functions.exponential(x, norm_noise, tau_noise, cost_noise), '', '', legend ='noise' )    
    n, bins = numpy.histogram(time_diff_up,  bins = n_bins, range = (0., 15.))
    bin_center = 0.5 * (bins[1:] + bins[:-1])    
    plot_functions.scatter_plot(bin_center, n, '', '' , fmt ='.b', legend = ' dati')    

    n = n - functions.exponential(bin_center, norm_noise, tau_noise, cost_noise)
    mask = (n > 0.)
    bin_center = bin_center[mask]
    n = n[mask]

    p0 = [10, 0.7 , 0.1, 2., 1.]
    plt.figure()
    plt.subplot(2,1,1)
    plt.title('dati sottratti')
    likelihood_fit(model = functions.two_expo_integral, n=n, bin_center = bin_center,
                   jacobian = functions.two_expo_integral_jacobian, n_bins = n_bins, x0 = p0,
                   fit_min = fit_min,  range_hist = (0., 15.) )

    p0 = [10, 0.7 , 0.1, 2., 1.]
    plt.figure()
    plt.title('dati non  sottratti')
    plt.subplot(2,1,1)
    likelihood_fit(model = functions.two_expo_integral, time_diff = time_diff_up,
                   jacobian = functions.two_expo_integral_jacobian, n_bins = n_bins, x0 = p0,
                   fit_min = 0.06,  range_hist = (0., 15.) )

    """
    plt.ion()
    plt.show()  
          

    
