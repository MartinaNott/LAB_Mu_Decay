import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import matplotlib.pyplot as plt
import numpy
import argparse 
import datetime
import time

import event_generator_functions
import geometry
import functions
import constants
import plot_functions

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('--number_events', '-n', default=1000, type=int, help='')
options_parser.add_argument('--output_file', '-o', default='None', type=str, help='')


if __name__ == '__main__' :   
    
    start_time = time.time()   
    options = vars(options_parser.parse_args())  
    N = options['number_events']
    output_file_events = options['output_file']

    #E_kin = event_generator_functions.energy_generator(N, constants.MUON_MASS, functions.spectrum, 1., 1.e5, -2.7, 300., 50.)        
    theta, phi = event_generator_functions.angle_generator(N, functions.dist_theta)   
    x_start, y_start = event_generator_functions.position_on_scintillator_generator(N, l = geometry.L, w = geometry.W)
    x_stop, y_stop, z = event_generator_functions.cosmic_rays_propagation(x_start, y_start, theta, phi, h = geometry.H, s_position = geometry.position)
    f = event_generator_functions.check_particle_position(x_stop, y_stop, l = geometry.L, w = geometry.W, s_position = geometry.position) 
    
    mask = numpy.full(N, True)


    data = numpy.vstack((x_start, y_start, theta, phi, x_stop, y_stop, f)).T    
    print("x_start, y_start, theta, phi, x_stop, x_stop, flag \n", data)                  
    alpha = numpy.sum(f)/N    
    alpha_err = numpy.sqrt(alpha * (1. - alpha)/N)
    alpha = alpha * 100
    alpha_err = alpha_err * 100
    print("Number of events hitting the scintillator/Total number of events on S1:", numpy.sum(f), "/", N, "=", alpha, '+-', alpha_err)

    if(output_file_events.endswith('.txt')): 
      header ='%s \nE[MeV], P [MeV], beta, x1[m], y1[m], theta, phi, x3[m], y3[m], flag\n' % datetime.datetime.now()
      fmt = ['%.4f', '%.4f', '%.4f', '%.4f', '%.6f', '%.2f', '%.2f', '%.4f', '%.6f', '%d']
      numpy.savetxt(output_file_events, numpy.transpose([E_kin[mask], P[mask], beta[mask], x_start[mask], y_start[mask], theta[mask], phi[mask], x_stop[mask], y_stop[mask], f[mask]]) , fmt=fmt, header=header)
      print("Output file saved!\n\n")   
    print("Time of execution: %s seconds " % (time.time() - start_time))

    plt.figure()    
    plt.subplot(2,1,1)
    plot_functions.plot_histogram(theta, 'theta', '', n_bins = 200, range = None)
    plt.subplot(2,1,2)
    plot_functions.plot_histogram(phi, 'phi', '', n_bins = 200, range = None)    

    range_hist = (0., 30.)           
    plt.figure()
    plt.subplot(2,2,1)
    plot_functions.plot_histogram(x_start, 'x_start', '', n_bins = 200, range = range_hist)          
    plt.subplot(2,2,2)      
    plot_functions.plot_histogram(y_start, 'y_start', '', n_bins = 200, range = range_hist)          
    plt.subplot(2,2,3)
    plot_functions.plot_histogram(x_stop[mask], 'x_stop', '', n_bins = 200, range = range_hist)          
    plt.subplot(2,2,4)      
    plot_functions.plot_histogram(y_stop[mask], 'y_stop', '', n_bins = 200, range = range_hist)          
            
      
      
    plt.ion()
    plt.show()  
      
      
