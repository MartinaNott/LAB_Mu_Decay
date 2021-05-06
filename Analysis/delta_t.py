import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import functions 
import utilities

def plot_channel_histogram(ch, time, channel_start, channel_stop, save_fig = False, range_hist = (0., 20.)):
    index, channel_diff, time_diff = utilities.mask_array(ch, time, channel_start, channel_stop)   
    plt.figure()
    title = 'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(channel_diff)
    plot_functions.plot_histogram(time_diff * 1.e6, "time [$\mu$s]", "", bins = 30, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = False) 
    if save_fig == True:
        figlabel = 'dt_%d_%d.pdf' % (channel_start, channel_stop)
        plt.savefig('%s' % figlabel , format = 'pdf')

    return channel_diff, time_diff 
    
   
description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', default='None', type=str, help='input_file')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file = options['input_file']
    
    ch, time = numpy.loadtxt(data_file, unpack=True)
    
    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 4)
    tau_short = 1/0.88
    tau_long = 1/2.2
    x = numpy.linspace(0., 20., 500)

    plt.plot(x, functions.exponential(x, 100, tau_long), '.r')    
    plt.plot(x, functions.two_expo(x, 100., tau_short, 100., tau_long), '.g')    
    channel_diff, time_diff = plot_channel_histogram(ch, time, 1, 3) 

    plt.plot(x, functions.exponential(x, 100, tau_long), '.r')






    
    plot_channel_histogram(ch, time, 4, 3)
    plot_channel_histogram(ch, time, 1, 1, range_hist = (0., 1.))
    plot_channel_histogram(ch, time, 4, 4)
    plot_channel_histogram(ch, time, 3, 3)
    plot_channel_histogram(ch, time, 3, 4)
    plot_channel_histogram(ch, time, 3, 1)
    plot_channel_histogram(ch, time, 4, 1)    
    
       

    
    
plt.ion()
plt.show()


