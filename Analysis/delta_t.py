import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import functions 
import utilities

def plot_channel_histogram(ch, time, channel_start, channel_stop, fit_function = None, param_names = None , param_units = None, p0 = None,  bounds = (-numpy.inf, numpy.inf), save_fig = False, range_hist = None, x_min = 0., label = ''):
    index, channel_diff, time_diff = utilities.mask_array(ch, time, channel_start, channel_stop)   
    time_diff = time_diff * 1.e6
    
    plt.figure()
    title = 'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(channel_diff)
    bins, n, dn = plot_functions.plot_histogram(time_diff, "time [$\mu$s]", "", n_bins = 30, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True) 
      
    if fit_function is not None:
        plot_functions.fit_histogram(bins, n, dn, param_names, param_units, fit_function = fit_function, p0 = p0, bounds = bounds, x_min = x_min)
  
    if save_fig == True:
        figlabel = 'dt_%d_%d_%s_piombo.pdf' % (channel_start, channel_stop, label)
        plt.savefig('%s' % figlabel , format = 'pdf')
    print("\n\n")
    return channel_diff, time_diff 
   
description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', default=None, type=str, help='input_file')
options_parser.add_argument('-save_fig', '-s', default=False, action='store_true', help='save fig')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file = options['input_file']
    save_fig = options['save_fig']
    ch, time = numpy.loadtxt(data_file, unpack=True)
    
    """    
    #fit 14 13 con doppia exp
    param_names = ['a_short', 'm_short', 'a_long', 'm_long', 'costant']
    param_units = ['$\mu ^-1$s', '$\mu$s', '$\mu ^-1$s', '$\mu$s', '$\mu ^-1$s']   
    
    p0 = [0.05, 0.88, 0.05, 2.2, 0.008]
    bounds =  (0.0, 0.1, 0., 2., 0.), (numpy.inf, 1.0, numpy.inf, 2.5 , 1.)
    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 4, fit_function = functions.two_expo, param_names = param_names, param_units = param_units, p0 = p0 , bounds = bounds, x_min = 0.5, range_hist = (0., 20.), save_fig=save_fig, label ='2exp')
       
    p0 = [0.05, 0.88, 0.05, 2.2, 0.008]
    bounds =  (0.0, 0.1, 0., 2., 0.), (numpy.inf, 1.0, numpy.inf, 2.5 , 1.)
    channel_diff, time_diff2 = plot_channel_histogram(ch, time, 1, 3, fit_function = functions.two_expo, param_names = param_names, param_units = param_units, p0 = p0, bounds = bounds, x_min = 0.2, range_hist = (0., 20.), save_fig=save_fig, label ='2exp') 
    
    #fit 14 13 con una sola exp
    param_names = ['a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '']    
    p0 = None
    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 4, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0 ,  x_min = 0.5, range_hist = (0., 20.), save_fig=save_fig)

    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 3, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0, x_min = 0.2, range_hist = (0., 20.),save_fig=save_fig) 

    p0 = [1 *1.e-15, 0.03 * 1.e6, 0.00001]    
    #plot_channel_histogram(ch, time, 1, 1, fit_function = functions.exponential, param_names = param_names, param_units = param_units, range_hist = (0., 1000000.), p0 = p0)
    
    plot_channel_histogram(ch, time, 4, 3, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 4, 4, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 3, 3, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 3, 4, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 3, 1 , range_hist = (0., 1000000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 4, 1, range_hist = (0., 3), save_fig=save_fig)    
    """   
 
     #fit 14 13 con doppia exp
    param_names = ['a_short', 'm_short', 'a_long', 'm_long', 'costant']
    param_units = ['$\mu ^-1$s', '$\mu$s', '$\mu ^-1$s', '$\mu$s', '$\mu ^-1$s']   
    
    p0 = [0.05, 0.88, 0.05, 2.2, 0.008]
    bounds =  (0.0, 0.1, 0., 2., 0.), (numpy.inf, 1.0, numpy.inf, 2.5 , 1.)
    channel_diff, time_diff = plot_channel_histogram(ch, time, 5, 7, fit_function = functions.two_expo, param_names = param_names, param_units = param_units, p0 = p0 , bounds = bounds, x_min = 0.5, range_hist = (0., 20.), save_fig=save_fig, label ='2exp')
       
    p0 = [0.05, 0.88, 0.05, 2.2, 0.008]
    bounds =  (0.0, 0.1, 0., 2., 0.), (numpy.inf, 1.0, numpy.inf, 2.5 , 1.)
    channel_diff, time_diff2 = plot_channel_histogram(ch, time, 5, 8, fit_function = functions.two_expo, param_names = param_names, param_units = param_units, p0 = p0, bounds = bounds, x_min = 0.2, range_hist = (0., 20.), save_fig=save_fig, label ='2exp') 
    
    #fit 14 13 con una sola exp
    param_names = ['a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '']    
    p0 = None
    channel_diff, time_diff = plot_channel_histogram(ch, time, 5, 8, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0 ,  x_min = 0.5, range_hist = (0., 20.), save_fig=save_fig)

    channel_diff, time_diff = plot_channel_histogram(ch, time, 5, 7, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0, x_min = 0.2, range_hist = (0., 20.),save_fig=save_fig) 

    p0 = [1 *1.e-15, 0.03 * 1.e6, 0.00001]    
    #plot_channel_histogram(ch, time, 1, 1, fit_function = functions.exponential, param_names = param_names, param_units = param_units, range_hist = (0., 1000000.), p0 = p0)
    
    plot_channel_histogram(ch, time, 8, 6, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 8, 8, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 7, 6, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 7, 8, range_hist = (0., 5000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 7, 5 , range_hist = (0., 1000000), save_fig=save_fig)
    plot_channel_histogram(ch, time, 8, 5, range_hist = (0., 3), save_fig=save_fig)    
         
         
    plt.ion()
    plt.show()


