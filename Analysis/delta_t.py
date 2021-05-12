import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import functions 
import utilities

def plot_channel_histogram(ch, time, channel_start, channel_stop, fit_function = None, param_names = None , param_units = None, p0 = None,  bounds = (-numpy.inf, numpy.inf), save_fig = False, range_hist = None, x_min = 0., label = ''):
    print('\n')
    index, channel_diff, time_diff = utilities.mask_array(ch, time, channel_start, channel_stop)   
    time_diff = time_diff * 1.e6
    
    title = 'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(channel_diff)
    bins, n, dn = plot_functions.plot_histogram(time_diff, "time [$\mu$s]", "", n_bins = 100, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True)      
    if fit_function is not None:
        opt, pcov = plot_functions.fit_histogram(bins, n, dn, param_names, param_units, fit_function = fit_function, p0 = p0, bounds = bounds, x_min = x_min)
        l_likelihood = utilities.log_likelihood(bins, n, dn, fit_function, *opt)
	#if save_fig == True:
	#    figlabel = 'dt_%d_%d_%s.pdf' % (channel_start, channel_stop, label)
        #    plt.savefig('%s' % figlabel , format = 'pdf')
        
        return channel_diff, time_diff , l_likelihood          
    return channel_diff, time_diff 
   
description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', default=None, type=str, help='input_file')
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
    ch, time = numpy.loadtxt(data_file, unpack=True)

    param_names_2exp = ['norm', 'fraction', 'm_short', 'm_long', 'costant']
    param_units_2exp = ['$\mu ^-1$s', '', '$\mu$s', '$\mu$s', '$\mu ^-1$s']  
    param_names = ['a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '']      
    
    p0 = [0.05, 0.5, 0.88, 2.2, 0.008]
    bounds =  (0.0, 0.0, 0.1, 1.5, 0.), (numpy.inf, 1., 1.5, 5., 1.)
    x_min = 0.0
    x_max = 20.
    
    plt.figure()        
    channel_diff_up, time_diff_up, l_likelihood_2exp = plot_channel_histogram(ch, time, ch_start, ch_stop_up, fit_function = functions.two_expo, param_names = param_names_2exp, param_units = param_units_2exp, p0 = p0 , bounds = bounds, x_min = x_min, range_hist = (0., 20.), save_fig=save_fig)       
    channel_diff, time_diff, l_likelihood_exp = plot_channel_histogram(ch, time, ch_start, ch_stop_up, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = None ,  x_min = x_min, range_hist = (0., 20.), save_fig=save_fig)       
    test = utilities.ll_ratio_test_stat(l_likelihood_2exp, l_likelihood_exp)
    print("test: ", test)
       
    plt.figure()        
    channel_diff_down, time_diff_down, l_likelihood_2exp = plot_channel_histogram(ch, time, ch_start, ch_stop_down, fit_function = functions.two_expo, param_names = param_names_2exp, param_units = param_units_2exp, p0 = p0, bounds = bounds, x_min = x_min, range_hist = (0., 20.), save_fig=save_fig)      
    channel_diff, time_diff, l_likelihood_exp = plot_channel_histogram(ch, time, ch_start, ch_stop_down, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = None, x_min = x_min, range_hist = (0., 20.), save_fig=save_fig) 

    test = utilities.ll_ratio_test_stat(l_likelihood_2exp, l_likelihood_exp)
    print("test: ", test)
    print("--------------------\n")

    ch_stop = numpy.concatenate((channel_diff_up, channel_diff_down)) 
    time_stop = numpy.concatenate((time_diff_up, time_diff_down)) 
    plt.figure()
    title = ''#'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(ch_stop)
    bins, n, dn = plot_functions.plot_histogram(time_stop, "time [$\mu$s]", "", n_bins = 120, range = (0., 20.), title = title, legend = legend, fmt = '.b', as_scatter = True) 
    plot_functions.fit_histogram(bins, n, dn, param_names_2exp, param_units_2exp, fit_function = functions.two_expo, p0 = p0, bounds = bounds, x_min = x_min, x_max = x_max) 
               

    plt.figure()    
    p0 = [1 *1.e-15, 0.03 * 1.e6, 0.00001]    
    plot_channel_histogram(ch, time, ch_start, ch_start, range_hist = (0., 500000.)) #fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0)

    plt.figure()    
    plt.subplot(2, 3, 1)
    plot_channel_histogram(ch, time, ch_stop_down, ch_stop_up, range_hist = (0., 5000), save_fig=save_fig)
    plt.subplot(2, 3, 2)
    plot_channel_histogram(ch, time, ch_stop_down, ch_stop_down, range_hist = (0., 5000), save_fig=save_fig)
    plt.subplot(2, 3, 3)
    plot_channel_histogram(ch, time, ch_stop_up, ch_stop_up, range_hist = (0., 5000), save_fig=save_fig)
    plt.subplot(2, 3, 4)
    plot_channel_histogram(ch, time, ch_stop_up, ch_stop_down, range_hist = (0., 5000), save_fig=save_fig)
    plt.subplot(2, 3, 5)
    plot_channel_histogram(ch, time, ch_stop_up, ch_start , range_hist = (0., 1000000), save_fig=save_fig)
    plt.subplot(2, 3, 6)
    plot_channel_histogram(ch, time, ch_stop_down, ch_start, range_hist = (0., 3), save_fig=save_fig)    


    plt.ion()
    plt.show()


