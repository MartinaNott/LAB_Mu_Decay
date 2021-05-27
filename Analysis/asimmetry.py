import sys
sys.path.insert(1, '/home/elo/ele/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse

import plot_functions
import functions
import utilities
import two_expo_fit

def calculate_asimmetry(time_diff_up, time_diff_down, n_bins, range_hist, x_min):
    mask = time_diff_up > x_min
    time_diff_up = time_diff_up[mask]
    mask = time_diff_down > x_min
    time_diff_down = time_diff_down[mask]  
    
    epsilon = len(time_diff_up)/len(time_diff_down)
    n_up, bins_up = numpy.histogram(time_diff_up,  bins = n_bins, range = range_hist) 
    n_down, bins_down = numpy.histogram(time_diff_down,  bins = n_bins, range = range_hist) 
    n_down = n_down * epsilon
    print("epsilon = ", epsilon)
    bins = bins_down
    bins_center = 0.5 * (bins[1:] + bins[:-1])    
    mask = (n_down + n_up) > 0.
    n_down = n_down[mask]
    n_up = n_up[mask]
    dn_up = numpy.sqrt(n_up)
    dn_down = numpy.sqrt(n_down)
    bins_center = bins_center[mask]
    asimmetry = (n_up - n_down)/(n_up + n_down)
    #asimmetry_err = ((1. + epsilon)/(n_up +n_down)**2) * numpy.sqrt(dn_up**2 + dn_down **2)
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

    x_min = 0.5
    x_max = 20.
    n_bins = 40
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch, time, ch_start, ch_stop_up)   
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch, time, ch_start, ch_stop_down)   
    range_hist = (time_diff_up[time_diff_up > 0.].min(), 20.)
    
    plt.figure()        
    two_expo_fit.plot_channel_histogram(time_diff_up, ch_start, ch_stop_up, n_bins = n_bins, range_hist = range_hist, save_fig=save_fig) 
    two_expo_fit.plot_channel_histogram(time_diff_down, ch_start, ch_stop_down, n_bins = n_bins, range_hist = range_hist, save_fig=save_fig) 

    #AGGREGANDO I DATI: SOPRA E SOTTO    
    ch_stop = numpy.concatenate((channel_diff_up, channel_diff_down)) 
    time_stop = numpy.concatenate((time_diff_up, time_diff_down)) 
    plt.figure()
    title = ''
    legend = '%d' % len(ch_stop)
    bins, n, dn = plot_functions.plot_histogram(time_stop, "time [$\mu$s]", "", n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True) 


    #ASIMMETRIA
    asimmetry, asimmetry_err, bins_center = calculate_asimmetry(time_diff_up, time_diff_down, n_bins, range_hist, x_min)
    plt.figure()
    plt.subplot(3,1,1)
    plot_functions.scatter_plot(bins_center, asimmetry, 'dt [$\mu$s]', 'Asimmetry ', dy = asimmetry_err, title = '')
    p0 = [0.1, 3., 0.0, 0.1]
    bounds = (0., 0., -numpy.inf, -numpy.inf), (0.3, numpy.inf, 2 * numpy.pi, +numpy.inf)
    param_names = ['Amplitude', '$\omega$', '$\phi$', 'costant']
    param_units = ['', 'MHz', 'rad', '']
    opt_wave, pcov_wave = plot_functions.do_fit(bins_center, asimmetry, asimmetry_err, param_names, param_units, functions.wave, p0 = p0, bounds = bounds , x_min = x_min, x_max = x_max)
    plt.subplot(3,1,2)    
    plot_functions.scatter_plot(bins_center, asimmetry, 'dt [$\mu$s]', 'Asimmetry ', dy = asimmetry_err, title = '')    
    p0 = [0.1, 3., 0.0, 0.1, 0.1]
    bounds = (0., 0., -numpy.inf, -numpy.inf, -numpy.inf), (0.3, numpy.inf, 2 * numpy.pi, +numpy.inf, numpy.inf )
    param_names = ['Amplitude', '$\omega$', '$\phi$', 'costant', '$\mu s^-1$']
    param_units = ['', 'MHz', 'rad', '', 'MHz']
    opt_increasing_wave, pcov_increasing_wave  = plot_functions.do_fit(bins_center, asimmetry, asimmetry_err, param_names, param_units, functions.increasing_wave, p0 = p0, bounds = bounds , x_min = x_min, x_max = x_max)
    plt.subplot(3,1,3)
    plot_functions.scatter_plot(bins_center, asimmetry, 'dt [$\mu$s]', 'Asimmetry ', dy = asimmetry_err, title = '')
    opt_line, pcov_line = plot_functions.do_fit(bins_center, asimmetry, asimmetry_err, param_names = ['m', 'costant'], param_units=['MHz', ''], fit_function=functions.line, p0 = None , x_min = x_min, x_max = x_max) 

    #figlabel = 'asimmetria_ferro_magnetizzato.pdf'
    #plt.savefig('%s' % figlabel , format = 'pdf')

    plt.figure()
    plt.subplot(2, 1, 1)
    plot_functions.scatter_plot(bins_center, asimmetry, 'dt [$\mu$s]', 'Asimmetry ', dy = asimmetry_err, title = '') 
    opt_costant, pcov_costant = plot_functions.do_fit(bins_center, asimmetry, asimmetry_err, param_names = ['costant'], param_units=[''], fit_function=functions.costant, p0 = None , x_min = x_min, x_max = x_max)        
    opt_line, pcov_line = plot_functions.do_fit(bins_center, asimmetry, asimmetry_err, param_names = ['m', 'costant'], param_units=['MHz', ''], fit_function=functions.line, p0 = None , x_min = x_min, x_max = x_max) 
    plt.subplot(2, 1, 2)
    residui = (asimmetry - functions.line(bins_center, *opt_line))/asimmetry_err
    plot_functions.scatter_plot(bins_center, residui, 'dt [$\mu$s]', 'Residui ', dy =asimmetry_err/asimmetry_err, title = '')



    #likelihood test: 
    mask = (bins_center > x_min)
    asimmetry = asimmetry[mask]
    asimmetry_err = asimmetry_err[mask]
    bins_center = bins_center[mask]
    l_likelihood_wave = functions.gauss_log_likelihood(bins_center, asimmetry, asimmetry_err, functions.wave, *opt_wave)    
    l_likelihood_increasing_wave = functions.gauss_log_likelihood(bins_center, asimmetry, asimmetry_err, functions.increasing_wave, *opt_increasing_wave)  
    l_likelihood_costant = functions.gauss_log_likelihood(bins_center, asimmetry, asimmetry_err, functions.costant, *opt_costant)  
    l_likelihood_line = functions.gauss_log_likelihood(bins_center, asimmetry, asimmetry_err, functions.line, *opt_line)          
    print('l_likelihood_wave', l_likelihood_wave)
    print('l_likelihood_increasing_wave', l_likelihood_increasing_wave)
    print('l_likelihood_costant', l_likelihood_costant)        
    print('l_likelihood_line', l_likelihood_line)        
    test_cost_wave = functions.ll_ratio_test_stat(l_likelihood_wave, l_likelihood_costant)
    test_line_increasing_wave = functions.ll_ratio_test_stat(l_likelihood_increasing_wave, l_likelihood_line)
    test_wave_increasing_wave = functions.ll_ratio_test_stat(l_likelihood_increasing_wave, l_likelihood_wave)
    print('\n\ntest_cost_wave', test_cost_wave, '\ntest_wave_increasing_wave', test_wave_increasing_wave, '\ntest_line_increasing_wave', test_line_increasing_wave)

    """
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
    plt.subplot(2,1,1) 
    two_expo_fit.plot_channel_histogram(time_diff, ch_start, ch_start, n_bins = 200, range_hist = (0., 2.5))#, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0)  
    plt.subplot(2,1,2)    
    two_expo_fit.plot_channel_histogram(time_diff, ch_start, ch_start, n_bins = 200, range_hist = (0., 100000))#, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0)    
    """
    plt.ion()
    plt.show()
        
        
