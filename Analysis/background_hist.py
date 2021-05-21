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
options_parser.add_argument('-input_file_up', '-u', default=None, type=str, help='input_file_up_circuit')
options_parser.add_argument('-input_file_down', '-d', default=None, type=str, help='input_file_down_circuit')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file_up = options['input_file_up']
    data_file_down = options['input_file_down']
    
    ch_up, time_up = numpy.loadtxt(data_file_up, unpack=True)
    ch_down, time_down = numpy.loadtxt(data_file_down, unpack=True)
        
    plt.figure()
    plt.subplot(2, 1, 1)
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch_up, time_up, 1, 1)   
    mask = time_diff_up>0.
    time_diff_up = time_diff_up[mask]
    print(min(time_diff_up))
    two_expo_fit.plot_channel_histogram(time_diff_up, 1, 1, n_bins = 50, range_hist = (0, 20.))

    plt.subplot(2, 1, 2)    
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch_down, time_down, 5, 5)   
    mask = time_diff_down>0.
    time_diff_down = time_diff_down[mask]
    print(min(time_diff_down))
    two_expo_fit.plot_channel_histogram(time_diff_down, 5, 5, n_bins = 50, range_hist = (0, 20.))

    plt.ion()
    plt.show()
