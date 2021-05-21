import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse

import plot_functions
import functions
import utilities
import two_expo_fit

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', nargs='*' , default=None, type=str, help='input_file')
options_parser.add_argument('-save_fig', '-s', default=False, action='store_true', help='save fig')
options_parser.add_argument('-ch_start', '-start', default=None, type=int, help='ch_start')
options_parser.add_argument('-ch_stop_up', '-up', default=None, type=int, help='ch_stop_up')
options_parser.add_argument('-ch_stop_down', '-down', default=None, type=int, help='ch_stop_down')

if __name__ == '__main__' :
    options = vars(options_parser.parse_args())
    data_file = options['input_file']
    save_fig = options['save_fig']
    ch_start = options['ch_start']
    ch_stop_up = options['ch_stop_up']
    ch_stop_down = options['ch_stop_down']

    data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file])
    ch = data[0, :]
    time = data[1, :]

    param_names_2exp = ['norm', 'fraction', 'm_short', 'm_long', 'costant']
    param_units_2exp = ['$\mu ^-1$s', '', '$\mu$s', '$\mu$s', '$\mu ^-1$s']

    #PARAMETRI INIZIALI DEL FIT E BOUNDARIES DEI PARAMETRI PER I FIT CON DOPPIA ESPONENZIALE
    p0 = [2.e+6, 0.5, 0.088, 2.2, 0.008]
    bounds =  (0., 0., 0.060, 1.5, 0.), (numpy.inf, 1., 0.090, 5., 1.)
    x_min = 0.0 #0.045# 0.64 
    x_max = 20.
    n_bins_up = 100
    n_bins_down = 100
    n_bins = 200

    integral = functions.do_expo_integral()
     
    #FIT DEI DATI CON SOPRA VERSO L'ALTO       
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch, time, ch_start, ch_stop_up)   
    range_hist = (time_diff_up[time_diff_up > 0.].min(), x_max)

    plt.figure()        
    plt.subplot(2, 1, 1)
    l_likelihood = two_expo_fit.plot_channel_histogram(time_diff_up, ch_start, ch_stop_up, n_bins = n_bins_up, fit_function = integral, param_names = param_names_2exp, param_units = param_units_2exp, p0 = p0 , bounds = bounds, x_min = x_min, range_hist = range_hist, save_fig=save_fig)       

    
    #FIT DEI DATI CON SOPRA VERSO IL BASSO   
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch, time, ch_start, ch_stop_down)   
    range_hist = (time_diff_down[time_diff_down > 0.].min(), x_max)
    plt.figure()
    plt.subplot(2, 1, )
    l_likelihood = two_expo_fit.plot_channel_histogram(time_diff_down, ch_start, ch_stop_down, n_bins = n_bins_up, fit_function = integral, param_names = param_names_2exp, param_units = param_units_2exp, p0 = p0 , bounds = bounds, x_min = x_min, range_hist = range_hist, save_fig=save_fig)  
     
    #AGGREGANDO I DATI: SOPRA E SOTTO    
    ch_stop = numpy.concatenate((channel_diff_up, channel_diff_down)) 
    time_stop = numpy.concatenate((time_diff_up, time_diff_down)) 
    plt.figure()
    title = ''
    legend = '%d' % len(time_stop)
    bins, n, dn = plot_functions.plot_histogram(time_stop, "time [$\mu$s]", "", n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True) 
    plot_functions.do_fit(bins, n, dn, param_names_2exp, param_units_2exp, fit_function = integral, p0 = p0, bounds = bounds, x_min = x_min, x_max = x_max) 
    
    plt.ion()
    plt.show()

