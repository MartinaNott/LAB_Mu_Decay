import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import utilities

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', default='None', type=str, help='input_file')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file = options['input_file']
    
    ch, time = numpy.loadtxt(data_file, unpack=True)

    ch_diff = numpy.ediff1d(ch)
    time_diff = numpy.ediff1d(time) * 1.e6

    index, channel_diff, time_diff = utilities.mask_array(ch, time, 2, 2)
    
    plt.figure()
    plot_functions.plot_histogram(time_diff * 1.e6, "time [$\mu$s]", "", bins = None, range = (0., 20.), title = '2, 2', legend = '', fmt = '.b', as_scatter = False)
    print(index)
    
    index, channel_diff, time_diff = utilities.mask_array(ch, time, 1, 2)
    plt.figure()
    plot_functions.plot_histogram(time_diff * 1.e6, "time [$\mu$s]", "", bins = None, range = (0., 20.), title = '1, 2', legend = '', fmt = '.b', as_scatter = False)    



plt.ion()
plt.show()


