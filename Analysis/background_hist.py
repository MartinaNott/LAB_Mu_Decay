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
        
    range_hist = (0., 15.)
    n_bins = 150
    param_names = ['N', '$\\tau$', 'costant']
    param_units = ['', '$\mu$s', '']
    
    mask_start = ch_down == 5
    a = numpy.ones(len(ch_down))
    a = a[mask_start]
    print('start down', len(a))

    mask_start = ch_up == 5
    b = numpy.ones(len(ch_up))
    b = b[mask_start]
    print('start up', len(b))

    title = '%d start' % (len(b))
    plt.figure()
    plt.subplot(2, 1, 1)
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch_up, time_up, 1, 4)   
    mask = time_diff_up>0.
    time_diff_up = time_diff_up[mask]
    print(min(time_diff_up))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_up, 1, 4, n_bins = n_bins, range_hist = range_hist, title = title, label = 'stop up, ')
    plot_functions.do_fit(bins, n, dn, param_names, param_units, fit_function = functions.exponential, p0 = None)
    plt.xlim(-1., 15)
    plt.ylim(-1., 25)
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch_up, time_up, 1, 3)       
    plt.subplot(2, 1, 2)
    mask = time_diff_up>0.
    time_diff_up = time_diff_up[mask]
    print(min(time_diff_up))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_up, 1, 3, n_bins = n_bins, range_hist = range_hist, title = '', label = 'stop down, ')
    #plot_functions.do_fit(bins, n, dn, param_names, param_units, fit_function = functions.exponential, p0 = None)
    plt.xlim(-1., 15)
    plt.ylim(-1., 25)

    title = '%d start' % (len(a))
    plt.figure()
    plt.subplot(2, 1, 1)    
    plt.ylim(-1., 25)
    plt.xlim(-1., 15)
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch_down, time_down, 5, 7)   
    mask = time_diff_down>0.
    time_diff_down = time_diff_down[mask]
    mask_background = time_diff_down < 0.2
    a = numpy.ones(len(time_diff_down))
    a = a[mask_background]
    print('dati con delta t tra 0 e 200ns:', len(a))
        
    print(min(time_diff_down))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_down, 5, 7, n_bins = n_bins, range_hist = range_hist, title = title, label = 'stop up, ')
    print(bins, n, len(bins), len(n))
    
    plot_functions.do_fit(bins, n, dn, param_names, param_units, fit_function = functions.exponential, p0 = None, x_min = 0.030)
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch_down, time_down, 5, 8)   
    
    plt.subplot(2,1,2)
    plt.ylim(-1., 25)
    plt.xlim(-1., 15)
    mask = time_diff_down>0.
    time_diff_down = time_diff_down[mask]
    print(min(time_diff_down))
    bins, n, dn = two_expo_fit.plot_channel_histogram(time_diff_down, 5, 8, n_bins = n_bins, range_hist = range_hist, title = '', label = 'stop down, ')
        
 
    #correction down: 16  6 eventi nel primo bin e nel secondo con 8522 start
    
    plt.ion()
    plt.show()
