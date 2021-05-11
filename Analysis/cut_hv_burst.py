import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse 

import plot_functions
import functions 
import utilities

def plot_channel_in_time(ch_mask): 
    plt.figure('events_sum_vs_row')
    rows = numpy.arange(len(ch_mask))
    row_bins = numpy.linspace(0, len(rows), 501)
    n, bins = numpy.histogram(rows[ch_mask], bins=row_bins)
    plt.hist(bins[0:-1],  weights = n, bins = bins)
    return 

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', default='None', type=str, help='input_file')
options_parser.add_argument('-boundaries', '-b', default='None', type=str, help='boundaries')
options_parser.add_argument('-output_file', '-o', default='None', type=str, help='output_file')

if __name__ == '__main__' :   
    options = vars(options_parser.parse_args())  
    data_file = options['input_file']
    boundaries = options['boundaries']
    output_file = options['output_file']
    ch, time = numpy.loadtxt(data_file, unpack=True)
    line_inf, line_sup = numpy.loadtxt(boundaries, unpack=True)   
         
    utilities.find_hv_bursts(ch)

    ch_mask = ch >= 3
    plot_channel_in_time(ch_mask)
     
    #mask = numpy.ediff1d(time) < -1
    #idx = numpy.where(mask)      
    #for i in idx[0]:
    #	time[i+1:] += time[i] - time[i+1]
    
    mask = numpy.ones(len(ch), dtype = bool)    
    for i in range (0, len(line_sup)-1): 
        mask = utilities.remove_hv_bursts(int(line_inf[i]), int(line_sup[i]), mask)    
        i = i + 1 
    
    print("number of line before cutting:" , len(ch))
    ch = ch[mask]
    time = time[mask]
    print("number of line after cutting:" , len(ch))
    
    if  output_file != 'None': 
       	header ='ch_FPGA time[s]\n' 
        fmt = ['%d', '%.8f']
        numpy.savetxt(output_file, numpy.transpose([ch, time]) , fmt=fmt, header=header)    
    
    
    plt.ion()
    plt.show()
    
      
