import sys
sys.path.insert(1, '/home/elo/ele/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import functions 
import utilities
import two_expo_fit

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file_up', '-u', nargs='*', default=None, type=str, help='input_file_up_circuit')
options_parser.add_argument('-input_file_down', '-d', nargs='*', default=None, type=str, help='input_file_down_circuit')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file_up = options['input_file_up']
    data_file_down = options['input_file_down']

    data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file_up])
    ch_up = data[0, :]
    time_up = data[1, :]    
    data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file_down])
    ch_down = data[0, :]
    time_down = data[1, :]
        
    range_hist = (0., 2.)
    n_bins = 100
    param_names = ['a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '']
    
    plt.figure()
    plt.subplot(2, 1, 1)
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch_up, time_up, 1, 4)   
    mask = time_diff_up>0.
    time_diff_up = time_diff_up[mask]
    print(min(time_diff_up))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_up, 1, 4, n_bins = n_bins, range_hist = range_hist)
    #plot_functions.do_fit(bins, n, dn, param_names, param_units, fit_function = functions.exponential, p0 = None)
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch_up, time_up, 1, 3)       
    """
    plt.subplot(2, 1, 2)
    mask = time_diff_up>0.
    time_diff_up = time_diff_up[mask]
    print(min(time_diff_up))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_up, 1, 3, n_bins = n_bins, range_hist = range_hist)
    #plot_functions.do_fit(bins, n, dn, param_names, param_units, fit_function = functions.exponential, p0 = None)
    """    
    plt.subplot(2, 1, 2)    
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch_down, time_down, 5, 7)   
    mask = time_diff_down>0.
    time_diff_down = time_diff_down[mask]
    print(min(time_diff_down))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_down, 5, 7, n_bins = n_bins, range_hist = range_hist)
    #plot_functions.do_fit(bins, n, dn, param_names, param_units, fit_function = functions.exponential, p0 = None)
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch_down, time_down, 5, 8)   
    """
    plt.subplot(2,1,2)
    mask = time_diff_down>0.
    time_diff_down = time_diff_down[mask]
    print(min(time_diff_down))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_down, 5, 8, n_bins = n_bins, range_hist = range_hist)
    """    
    
    plt.ion()
    plt.show()
