import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt
import argparse

import plot_functions
import functions
import utilities

def plot_channel_histogram(time_diff, channel_start, channel_stop, n_bins, fit_function = None, param_names = None , param_units = None, p0 = None,  bounds = (-numpy.inf, numpy.inf), save_fig = False, range_hist = None, x_min = 0., label = '', ex_int = (numpy.inf, -numpy.inf)):
    print('\n')
    title = 'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(time_diff)
    bins, n, dn = plot_functions.plot_histogram(time_diff, "time [$\mu$s]", "", n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True)      
    if fit_function is not None:
        opt, pcov = plot_functions.fit_histogram(bins, n, dn, param_names, param_units, fit_function = fit_function, p0 = p0, bounds = bounds, x_min = x_min, ex_int = ex_int)
        l_likelihood = utilities.log_likelihood(bins, n, dn, fit_function, *opt)
	#if save_fig == True:
	#    figlabel = 'dt_%d_%d_%s.pdf' % (channel_start, channel_stop, label)
        #    plt.savefig('%s' % figlabel , format = 'pdf')
        
        return  l_likelihood          
    return 

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

    param_names_2exp = ['norm', 'fraction', 'm_short', 'm_long', 'costant']
    param_units_2exp = ['$\mu ^-1$s', '', '$\mu$s', '$\mu$s', '$\mu ^-1$s']
    param_names_2expo_gauss = ['norm', 'fraction', 'm_short', 'm_long', 'costant', 'norm', 'mean', 'sigma']
    param_units_2expo_gauss = ['$\mu ^-1$s', '', '$\mu$s', '$\mu$s', '$\mu ^-1$s', '', 'mus', 'mus']
    param_names = ['a_long', 'm_long', 'costant']
    param_units = ['1/$\mu$s', '$\mu$s', '']

    #PARAMETRI INIZIALI DEL FIT E BOUNDARIES DEI PARAMETRI PER I FIT CON DOPPIA ESPONENZIALE
    p0 = [1., 0.5, 0.088, 2.2, 0.008]
    bounds =  (0.0, 0.01, 0.02, 1.5, 0.), (numpy.inf, 0.999, 1.3, 5., 1.)
    x_min = 0.0 #0.045# 0.64 
    x_max = 20.
    n_bins_up = 80
    n_bins_down = 80
    n_bins = 80
       
    index, channel_diff_up, time_diff_up = utilities.mask_array(ch, time, ch_start, ch_stop_up)   
    range_hist = (time_diff_up[time_diff_up > 0.].min(), x_max)

    #FIT DEI DATI CON SOPRA VERSO L'ALTO
    plt.figure()        
    l_likelihood_2exp = plot_channel_histogram(time_diff_up, ch_start, ch_stop_up, n_bins = n_bins_up, fit_function = functions.two_expo, param_names = param_names_2exp, param_units = param_units_2exp, p0 = p0 , bounds = bounds, x_min = x_min, range_hist = range_hist, save_fig=save_fig)       
    l_likelihood_exp = plot_channel_histogram(time_diff_up, ch_start, ch_stop_up, n_bins = n_bins_up, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = None ,  x_min = x_min, range_hist = range_hist, save_fig=save_fig, ex_int=(1.4, 4.5))       
    test = utilities.ll_ratio_test_stat(l_likelihood_2exp, l_likelihood_exp)
    print("test: ", test)
    index, channel_diff_down, time_diff_down = utilities.mask_array(ch, time, ch_start, ch_stop_down)   
    range_hist = (time_diff_down[time_diff_down > 0.].min(), x_max)


    #FIT DEI DATI CON SOPRA VERSO IL BASSO   
    x_min = 0. #0.045# 0.64 
    plt.figure()        
    l_likelihood_2exp = plot_channel_histogram(time_diff_down, ch_start, ch_stop_down, n_bins = n_bins_down, fit_function = functions.two_expo, param_names = param_names_2exp, param_units = param_units_2exp, p0 = p0, bounds = bounds, x_min = x_min, range_hist = range_hist, save_fig=save_fig, ex_int=(1.4, 4.5))      
    l_likelihood_exp = plot_channel_histogram(time_diff_down, ch_start, ch_stop_down, n_bins = n_bins_down, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = None, x_min = x_min, range_hist = range_hist, save_fig=save_fig) 


    test = utilities.ll_ratio_test_stat(l_likelihood_2exp, l_likelihood_exp)
    print("test: ", test)
    print("-------\n")

    #AGGREGANDO I DATI: SOPRA E SOTTO    
    ch_stop = numpy.concatenate((channel_diff_up, channel_diff_down)) 
    time_stop = numpy.concatenate((time_diff_up, time_diff_down)) 
    x_min = 0. #0.045# 0.64     
    plt.figure()
    title = ''#'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(ch_stop)
    bins, n, dn = plot_functions.plot_histogram(time_stop, "time [$\mu$s]", "", n_bins = n_bins, range = range_hist, title = title, legend = legend, fmt = '.b', as_scatter = True) 
    plot_functions.fit_histogram(bins, n, dn, param_names_2exp, param_units_2exp, fit_function = functions.two_expo, p0 = p0, bounds = bounds, x_min = x_min, x_max = x_max, ex_int=(1.4, 4.5)) 
    #plot_functions.fit_histogram(bins, n, dn, param_names = param_names, param_units = param_units, fit_function = functions.exponential, p0 = None, bounds = (-numpy.inf, numpy.inf), x_min = x_min, x_max = x_max)
    
    asimmetry, asimmetry_err, bins = calculate_asimmetry(time_diff_up, time_diff_down, n_bins, range_hist)
    plot_functions.scatter_plot(bins, asimmetry, 'dt', 'A', dy = asimmetry_err, title = '')
    
    
    #ALCUNI PLOT DI TEST PER VEDERE CORRELAZIONI, AFTERPULSE, E COSE CHE SUCCEDONO
    plt.figure()
    n_bins = 20
    plt.subplot(2, 3, 1)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_down, ch_stop_up)   
    plot_channel_histogram(time_diff, ch_stop_down, ch_stop_up, n_bins = n_bins, range_hist = (-0.1, 0.1), save_fig=save_fig)
    plt.subplot(2, 3, 2)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_down, ch_stop_down)   
    plot_channel_histogram(time_diff, ch_stop_down, ch_stop_down, n_bins = n_bins, save_fig=save_fig)
    plt.subplot(2, 3, 3)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_up, ch_stop_up)   
    plot_channel_histogram(time_diff, ch_stop_up, ch_stop_up, n_bins = n_bins, save_fig=save_fig)
    plt.subplot(2, 3, 4)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_up, ch_stop_down)   
    plot_channel_histogram(time_diff, ch_stop_up, ch_stop_down, n_bins = n_bins, save_fig=save_fig)
    plt.subplot(2, 3, 5)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_up, ch_start)   
    plot_channel_histogram(time_diff, ch_stop_up, ch_start, n_bins = n_bins, range_hist = (0., 10000), save_fig=save_fig)
    plt.subplot(2, 3, 6)
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_stop_down, ch_start)   
    plot_channel_histogram(time_diff, ch_stop_down, ch_start, n_bins = n_bins, save_fig=save_fig)    

    """          
    #DELTA T RELATIVO AI SEGNALI DI START (PER SAPERE PIU O MENO IL RATE DEI MUONI CHE SI FERMANO)
    p0 = [1000., 1.e5, 0.008]
    index, channel_diff, time_diff = utilities.mask_array(ch, time, ch_start, ch_start)   
    plt.figure()    
    plot_channel_histogram(time_diff, ch_start, ch_start, n_bins = n_bins, range_hist = (0., 1.))#, fit_function = functions.exponential, param_names = param_names, param_units = param_units, p0 = p0)    
    """    
    """
    #PLOT DI DUE ESPONENZIALI PIU UNA GAUSSIANA PICCOLINA (PER MODELLIZZARE EVENTUALMENTE IL BUMP  CHE SI VEDE IN ALCUNI FILE A 2.2 MICROSECONDI)
    p0 = [0.05, 0.5, 0.08, 2.2, 0.008, 0.05, 2.8, 0.2]
    bounds =  (0.0, 0.0, 0.05, 1.5, 0., 0.0, 0., 0.), (numpy.inf, 1., 1.2, 5., 1., numpy.inf, numpy.inf, numpy.inf)
    plt.figure()
    title = ''#'start:%d, stop:%d' %(channel_start, channel_stop)
    legend = '%d' % len(ch_stop)
    bins, n, dn = plot_functions.plot_histogram(time_stop, "time [$\mu$s]", "", n_bins = 100, range = (0., 20.), title = title, legend = legend, fmt = '.b', as_scatter = True) 
    plot_functions.fit_histogram(bins, n, dn, param_names_2expo_gauss, param_units_2expo_gauss, fit_function = functions.two_expo_gauss, p0 = p0, bounds = bounds, x_min = x_min, x_max = x_max) 
    """  
    

    plt.ion()
    plt.show()
