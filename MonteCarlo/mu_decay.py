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



def likelihood_fit(model, param_names, param_units, time_diff=None, bin_center=None, n=None, jacobian = None,  x0=None, bounds=None, title='', legend='', fit_min = -numpy.inf, range_hist = (0., 20), n_bins = 100, output_file = ''):
    
    plt.subplot(3,1,(1,2))
    if time_diff is not None: 
        bin_width = 100 #int(1000 * (range_hist[1]- range_hist[0])/n_bins)
        ylabel = 'ev./%d ns ' % (bin_width)
        bin_center, n, dn = plot_functions.plot_histogram(time_diff, "$\Delta t$ [$\mu$s]", ylabel, n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True)  
    else: 
        bin_width = int(1000 * (bin_center[2]- bin_center[1]))
        ylabel = 'ev./%d ns' % (bin_width)
        dn = numpy.sqrt(n)
        plot_functions.scatter_plot(bin_center, n, "$\Delta t$ [$\mu$s]", ylabel, dy = dn , title = title, fmt = '.b')    
    
    mask = bin_center >  fit_min   
    bin_center = bin_center[mask]
    n = n[mask]
    dn = dn[mask]
    minus_two_ll = functions.poisson_log_likelihood(bin_center, n, model)
    if jacobian is not None:
        jac = jacobian(bin_center, n)
    else:
        jac = None
    res = minimize(minus_two_ll, x0 = x0, bounds = bounds , method = 'BFGS', jac = jac, options={'disp' : True})
    opt = res.x 
    pcov = res.hess_inv
    L = minus_two_ll(opt)   
    param_results = plot_functions.fit_legend(opt, numpy.sqrt(pcov.diagonal()), param_names, param_units)
 
    plot_functions.scatter_plot(bin_center, model(bin_center, *opt), "$\Delta t$ [$\mu$s]", ylabel, fmt='-', legend = param_results, title = title)   
    
    plt.subplot(3,1,3)
    residuals = n - model(bin_center, *opt)
    plot_functions.scatter_plot(bin_center, residuals, "$\Delta t$ [$\mu$s]", "res", fmt='.')    
 
    if output_file is not '': 
        param_results = param_results + '\n\n'
        with open(output_file, 'a') as of:
            of.write(param_results)
    return L, bin_center, n, dn

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('--number_events', '-n', default=1000, type=int, help='')
options_parser.add_argument('--output_file', '-o', default='', type=str, help='')


if __name__ == '__main__' :   
    numpy.random.seed(7)
        
    start_time = time.time()   
    options = vars(options_parser.parse_args())  
    N = options['number_events']
    output_file = options['output_file']

    if output_file is not '':
        notes = 'Prima i risultati del fit del MC con noise, poi i risultati del fit dopo aver sottratto il noise.'
        with open(output_file, 'a') as of:
            of.write(notes)    

    param_names_2exp = ['norm', 'fraction', '$\\tau_s$', '$\\tau_l$', 'costant']
    param_units_2exp = ['$\mu s^-1$', '', '$\mu$s', '$\mu$s', '$\mu s^-1$']
    param_names_exp = ['N', '$\\tau$', 'costant']
    param_units_exp = ['$\mu s^-1$', '$\mu$s', '$\mu s^-1$']
    
    #norm, fraction, m_short, m_long, costant
    p0_gen = [1.586, 0.50054369 , 0.230, 2.17491405, 0.007, 0.06, 2.71, 0.4, 0.01, 4.5, 0.2, 0.021, 8., 0.2]
    n_bins = 150
    fit_min = 0.000


    norm_noise = 2.39 * 1366595/85227 #3291430/85227
    tau_noise = 0.106
    cost_noise = 1.2
    N_noise = int(100 * 1366595/85227)
    N_events = N - N_noise
    p0_noise_gen = [norm_noise, tau_noise, cost_noise]
    dt = mu_decay(N_events, functions.two_expo_n_gauss, 0., 100., *p0_gen)
    dt_noise = mu_decay(N_noise, functions.exponential, 0., 100., *p0_noise_gen)
    dt = numpy.concatenate((dt, dt_noise))


    #MONTE CARLO CON NOISE    
    title = 'Alluminio - %d eventi' %N
    plt.figure(figsize = [6.4, 7.])
    p0 = [1.44, 0.58054369 , 0.200, 2.197197491405, 32.]
    minus_two_ll, bin_center, n, dn = likelihood_fit( model = functions.two_expo_integral, param_names = param_names_2exp, param_units = param_units_2exp, 
                    time_diff= dt, jacobian = functions.two_expo_integral_jacobian, 
                    n_bins = n_bins, x0 = p0, bounds = None, fit_min = fit_min,  
                    range_hist = (0., 15.), title =title, output_file = output_file)

    plt.subplot(3,1,(1,2))
    bounds =  (0.0, 0.01, 0.02, 1.5, 0.), (numpy.inf, 1., 1.2, 5., 1000)
    opt, pcov = plot_functions.do_fit(bin_center, n, dn, param_names=param_names_2exp, param_units=param_units_2exp, fit_function = functions.two_expo_integral, p0 = p0, show=True, bounds = bounds, draw_on_points=True, output_file = output_file)

    #MONTE CARLO SOTTRAENDO NOISE        
    plt.figure(figsize = [6.4, 7.])
    n = n - functions.exponential(bin_center, norm_noise, tau_noise, cost_noise)
    mask = (n > 0.)
    bin_center = bin_center[mask]
    n = n[mask]
    dn = numpy.sqrt(n)
    minus_two_ll_two_expo, bin_center, n, dn = likelihood_fit( model = functions.two_expo_integral, param_names = param_names_2exp, param_units = param_units_2exp, 
                    bin_center=bin_center, n = n, jacobian = 
                    functions.two_expo_integral_jacobian, 
                    n_bins = n_bins, x0 = p0, bounds = None, fit_min = fit_min,  
                    range_hist = (0., 15.), title =title, output_file = output_file)
    
    plt.subplot(3,1,(1,2))
    opt, pcov = plot_functions.do_fit(bin_center, n, dn, param_names=param_names_2exp, param_units=param_units_2exp, fit_function = functions.two_expo_integral, p0 = p0, show=True, bounds = bounds, draw_on_points=True, output_file = output_file)

    #MONTE CARLO SOTTRAENDO NOISE fittato con un exp   
    plt.figure(figsize = [6.4, 7.])
    p0_one_exp = [6000., 2.2, 0.7]
    minus_two_ll_expo, bin_center, n, dn = likelihood_fit( model = functions.expo_integral, param_names=param_names_exp, param_units=param_units_exp,
                    bin_center=bin_center, n = n, jacobian = 
                    functions.expo_integral_jacobian, 
                    n_bins = n_bins, x0 = p0_one_exp, bounds = None, fit_min = fit_min,  
                    range_hist = (0., 15.), title =title, output_file = output_file)
    
    plt.subplot(3,1,(1,2))
    opt, pcov = plot_functions.do_fit(bin_center, n, dn, param_names=param_names_exp, param_units=param_units_exp, fit_function = functions.expo_integral, p0 = p0_one_exp, show=True, draw_on_points=True, output_file = output_file)

    print('\n\n\nminus_two_ll_two_expo, minus_two_ll_expo', minus_two_ll_two_expo, minus_two_ll_expo)
    test = minus_two_ll_expo - minus_two_ll_two_expo
    print("\ntest MC dopo aver sottratto il noise: ", test, '\n')

    
    
    data_file = 'data/run15_mudecay.dat'
    ch, time = numpy.loadtxt(data_file, unpack=True)
    #data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file])
    #ch = data[0, :]
    #time = data[1, :]
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch, time, 5, 7)   
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch, time, 5, 8)   

    mask_start = ch == 5
    a = numpy.ones(len(ch))
    a = a[mask_start]
    print('START down', len(a)) 
    print('events ', len(time_diff_up)+ len(time_diff_down))

    data_file = 'data/run11_mudecay.dat'
    ch, time = numpy.loadtxt(data_file, unpack=True)
    mask_start = ch == 1
    a = numpy.ones(len(ch))
    a = a[mask_start]
    print('START down', len(a)) 
    index, channel_diff_up_11, time_diff_up_11 = utilities.mask_array(ch, time, 1, 4)   
    index, channel_diff_down_11, time_diff_down_11 = utilities.mask_array(ch, time, 1, 3)   
    time_diff_up = numpy.concatenate((time_diff_up, time_diff_up_11))
    time_diff_down = numpy.concatenate((time_diff_down, time_diff_down_11))    
    print('events ', len(time_diff_up)+ len(time_diff_down))
    time_diff = numpy.concatenate((time_diff_up, time_diff_down))

 
    #DATI 
    title = 'Alluminio - %d eventi' %len(time_diff_up)
    plt.figure(figsize = [6.4, 7.])
    plt.subplot(3,1,(1,2))
    minus_two_ll, bin_center, n, dn = likelihood_fit(model = functions.two_expo_integral, param_names = param_names_2exp, param_units = param_units_2exp, 
                   time_diff = time_diff_up,
                   jacobian = functions.two_expo_integral_jacobian, n_bins = n_bins, x0 = p0,
                   fit_min = fit_min,  range_hist = (0., 15.), output_file = output_file , title=title)
    plt.subplot(3,1,(1,2))
    opt, pcov = plot_functions.do_fit(bin_center, n, dn, param_names=param_names_2exp, param_units=param_units_2exp, fit_function = functions.two_expo_integral, p0 = p0, show=True, bounds = bounds, draw_on_points=True, output_file = output_file)

    n = n - functions.exponential(bin_center, norm_noise, tau_noise, cost_noise)
    mask = (n > 0.)
    bin_center = bin_center[mask]
    n = n[mask]
    
    #DATI SOTTRATTI DAL NOISE
    plt.figure(figsize = [6.4, 7.])
    plt.subplot(3,1,(1,2))
    minus_two_ll, bin_center, n, dn = likelihood_fit(model = functions.two_expo_integral, param_names = param_names_2exp, param_units = param_units_2exp,  n=n, bin_center = bin_center,
                   jacobian = functions.two_expo_integral_jacobian, n_bins = n_bins, x0 = p0,
                   fit_min = fit_min,  range_hist = (0., 15.), output_file = output_file , title = title)
    print('len(n), len(dn)', len(n), len(dn))
    plt.subplot(3,1,(1,2))
    opt, pcov = plot_functions.do_fit(bin_center, n, dn, param_names=param_names_2exp, param_units=param_units_2exp, fit_function = functions.two_expo_integral, p0 = p0, show=True, bounds = bounds, draw_on_points=True, output_file = output_file, x_min = fit_min)
    
    
    plt.ion()
    plt.show()  
          

    
