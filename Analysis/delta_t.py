import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import functions 
import utilities

def plot_channel_histogram(ch, time, channel_start, channel_stop, fit_function, param_names, param_units, p0 = None,  bounds = (-numpy.inf, numpy.inf), save_fig = False, range_hist = (0., 20.), x_min = 0.):
    index, channel_diff, time_diff = utilities.mask_array(ch, time, channel_start, channel_stop)   
    time_diff = time_diff * 1.e6
    
    plt.figure()
    title = 'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(channel_diff)
    bins, n, dn = plot_functions.plot_histogram(time_diff, "time [$\mu$s]", "", n_bins = 120, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True) 
      
    plot_functions.fit_histogram(bins, n, dn, param_names, param_units, fit_function = fit_function, p0 = p0, bounds = bounds, x_min = x_min)
  
    if save_fig == True:
        figlabel = 'dt_%d_%d.pdf' % (channel_start, channel_stop)
        plt.savefig('%s' % figlabel , format = 'pdf')
    print("\n\n")
    return channel_diff, time_diff 

   
description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', default='None', type=str, help='input_file')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file = options['input_file']
    
    ch, time = numpy.loadtxt(data_file, unpack=True)
    
    param_names = ['a_short', 'm_short', 'a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '1/$\mu$s', '$\mu$s', '']
    
    p0 = [0.05, 0.88, 0.02, 2.2, 0.008]
    bounds =  (0.045, 0.75, 0., 2., 0.), (0.06, 1.0, 1., 2.5 , 1.)
    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 4, fit_function = functions.two_expo, param_names = param_names, param_units = param_units, p0 = p0 , bounds = bounds, x_min = 0.5)
    
    
    
    p0 = [0.05, 0.88, 0.02, 2.2, 0.008]
    bounds =  (0., 0.75, 0., 1.52, 0.), (1., 1.5, 1., 2.5 , 0.009)
    channel_diff, time_diff2 = plot_channel_histogram(ch, time, 1, 3, fit_function = functions.two_expo, param_names = param_names, param_units = param_units, p0 = p0, bounds = bounds, x_min = 0.2) 
    
    
    time_diff = numpy.concatenate(time_diff,time_diff2)
    plt.figure()
    title = 'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(channel_diff)
    bins, n, dn = plot_functions.plot_histogram(time_diff, "time [$\mu$s]", "", n_bins = 120, range = (0., 20.), title = title, legend = legend, fmt = '.b', as_scatter = True) 
      
    plot_functions.fit_histogram(bins, n, dn, param_names, param_units, fit_function = two_expo, p0 = p0, bounds = bounds, x_min = 0.2)    
    


    param_names = ['a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '']    
    p0 = None
    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 4, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0 ,  x_min = 0.5)

    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 3, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0, x_min = 0.2) 

    p0 = [1 *1.e-15, 0.03 * 1.e6, 0.00001]    
    #plot_channel_histogram(ch, time, 1, 1, fit_function = functions.exponential, param_names = param_names, param_units = param_units, range_hist = (0., 1000000.), p0 = p0)
    
    """
    plot_channel_histogram(ch, time, 4, 3, range_hist = (0., 0.5))       
    plot_channel_histogram(ch, time, 4, 4, range_hist = (0., 10.))
    plot_channel_histogram(ch, time, 3, 3)
    plot_channel_histogram(ch, time, 3, 4)
    plot_channel_histogram(ch, time, 3, 1)
    plot_channel_histogram(ch, time, 4, 1)    
    """
       

    
    
plt.ion()
plt.show()


