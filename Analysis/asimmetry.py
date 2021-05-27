import sys
sys.path.insert(1, '/home/elo/ele/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse

import plot_functions
import functions
import utilities
import two_expo_fit

def calculate_asimmetry(time_diff_up, time_diff_down, n_bins, range_hist):
    n_up, bins_up = numpy.histogram(time_diff_up,  bins = n_bins, range = range_hist) 
    n_down, bins_down = numpy.histogram(time_diff_down,  bins = n_bins, range = range_hist) 
    bins = bins_down
    bins_center = 0.5 * (bins[1:] + bins[:-1])    
    mask = (n_down + n_up) > 0.
    n_down = n_down[mask]
    n_up = n_up[mask]
    dn_up = numpy.sqrt(n_up)
    dn_down = numpy.sqrt(n_down)
    bins_center = bins_center[mask]
    asimmetry = (n_up - n_down)/(n_up + n_down)
    asimmetry_err = (2 / (n_up +n_down)**2) * numpy.sqrt((n_down * dn_up)**2 + (n_up * dn_down)**2) 
    return asimmetry, asimmetry_err,  bins_center
    

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('-input_file', '-f', nargs='*' , default=None, type=str, help='input_file')
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

    data = numpy.hstack([numpy.loadtxt(_file, unpack=True) for _file in data_file])
    ch = data[0, :]
    time = data[1, :]

    x_min = 0.3  
    x_max = 8.
    n_bins = 80
       
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch, time, ch_start, ch_stop_up)   
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch, time, ch_start, ch_stop_down)   
    range_hist = (time_diff_up[time_diff_up > 0.].min(), 10.)
    range_hist = (time_diff_down[time_diff_down > 0.].min(), 10.)
    
    plt.figure()        
    plt.subplot(2, 1, 1)
    two_expo_fit.plot_channel_histogram(time_diff_up, ch_start, ch_stop_up, n_bins = n_bins, range_hist = range_hist, save_fig=save_fig) 
    plt.subplot(2, 1, 2)
    two_expo_fit.plot_channel_histogram(time_diff_down, ch_start, ch_stop_down, n_bins = n_bins, range_hist = range_hist, save_fig=save_fig) 

    #AGGREGANDO I DATI: SOPRA E SOTTO    
    ch_stop = numpy.concatenate((channel_diff_up, channel_diff_down)) 
    time_stop = numpy.concatenate((time_diff_up, time_diff_down)) 
    plt.figure()
    title = ''
    legend = '%d' % len(ch_stop)
    bins, n, dn = plot_functions.plot_histogram(time_stop, "time [$\mu$s]", "", n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True) 

    #ASIMMETRIA
    asimmetry, asimmetry_err, bins_center = calculate_asimmetry(time_diff_up, time_diff_down, n_bins, range_hist)
    plt.figure()
    plot_functions.scatter_plot(bins_center, asimmetry, 'dt [$\mu$s]', 'Asimmetry ', dy = asimmetry_err, title = '')
    p0 = [0.1, 3., 0.0, 0.1]
    bounds = (0., 0., -numpy.inf, 0.), (0.3, numpy.inf, 2 * numpy.pi, 1. )
    param_names = ['Amplitude', '$\omega$', '$\phi$', '']
    param_units = ['', 'MHz', 'rad', '']
    #plot_functions.do_fit(bins_center, asimmetry, asimmetry_err, param_names, param_units, functions.wave, p0 = p0, bounds = bounds , x_min = x_min, x_max = x_max)
    
    
        
    #ALCUNI PLOT DI MONITORAGGIO
    plt.figure()
    n_bins = 20
    plt.subplot(2, 3, 1)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_down, ch_stop_up)   
    two_expo_fit.plot_channel_histogram(time_diff, ch_stop_down, ch_stop_up, n_bins = n_bins, range_hist = (-0.1, 0.1), save_fig=save_fig)
    plt.subplot(2, 3, 2)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_down, ch_stop_down)   
    two_expo_fit.plot_channel_histogram(time_diff, ch_stop_down, ch_stop_down, n_bins = n_bins, save_fig=save_fig)
    plt.subplot(2, 3, 3)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_up, ch_stop_up)   
    two_expo_fit.plot_channel_histogram(time_diff, ch_stop_up, ch_stop_up, n_bins = n_bins, save_fig=save_fig)
    plt.subplot(2, 3, 4)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_up, ch_stop_down)   
    two_expo_fit.plot_channel_histogram(time_diff, ch_stop_up, ch_stop_down, n_bins = n_bins, save_fig=save_fig)
    plt.subplot(2, 3, 5)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_up, ch_start)   
    two_expo_fit.plot_channel_histogram(time_diff, ch_stop_up, ch_start, n_bins = n_bins, range_hist = (0., 10000), save_fig=save_fig)
    plt.subplot(2, 3, 6)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_down, ch_start)   
    two_expo_fit.plot_channel_histogram(time_diff, ch_stop_down, ch_start, n_bins = n_bins, save_fig=save_fig)           
        
    #DELTA T RELATIVO AI SEGNALI DI START (PER SAPERE PIU O MENO IL RATE DEI MUONI CHE SI FERMANO)
    p0 = [1000., 1.e5, 0.008]
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_start, ch_start)   
    plt.figure()    
    two_expo_fit.plot_channel_histogram(time_diff, ch_start, ch_start, n_bins = n_bins, range_hist = (0., 1.))#, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0)  
    
    plt.ion()
    plt.show()
        
        
