import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse

import plot_functions
import functions
import utilities
import delta_t


description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', nargs='*' , default=None, type=str, help='input_file')
options_parser.add_argument('-save_fig', '-s', default=False, action='store_true', help='save fig')
options_parser.add_argument('-ch_start', '-start', default=None, type=int, help='ch_start')
options_parser.add_argument('-ch_middle', '-middle', default=None, type=int, help='')
options_parser.add_argument('-ch_stop', '-stop', default=None, type=int, help='')

if __name__ == '__main__' :
    options = vars(options_parser.parse_args())
    data_file = options['input_file']
    save_fig = options['save_fig']
    ch_start = options['ch_start']
    ch_middle = options['ch_middle']
    ch_stop = options['ch_stop']

    data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file])
    ch = data[0, :]
    time = data[1, :]


    index, channel_diff, time_diff = utilities.find_sequences_in_array(ch, time, ch_start, ch_middle, ch_stop)
    range_hist = (time_diff[time_diff > 0.].min(), 20.)
    
    plt.figure()        
    n_bins = 100
    delta_t.plot_channel_histogram(time_diff, ch_middle, ch_stop, n_bins = n_bins, range_hist = range_hist)  

    plt.ion()
    plt.show()
