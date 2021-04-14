import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
import argparse 
import datetime
import time

import event_generator_functions
import geometry
import functions
import constants

description = ''
options_parser = argparse.ArgumentParser(description = description)
options_parser.add_argument('--number_events', '-n', default=1000, type=int, help='')
options_parser.add_argument('--output_file', '-o', default='None', type=str, help='')


if __name__ == '__main__' :   
    
    start_time = time.time()   
    options = vars(options_parser.parse_args())  
    N = options['number_events']
    output_file_events = options['output_file']

    E_kin = event_generator_functions.energy_generator(N, constants.MUON_MASS, functions.spectrum, 1., 1.e5, -2.7, 300., 50.)        
    theta, phi = event_generator_functions.angle_generator(N, functions.dist_theta)   
    x_start, y_start = event_generator_functions.position_on_scintillator_generator(N, l = geometry.L, w = geometry.W)
    x_stop, y_stop, z = event_generator_functions.cosmic_rays_propagation(x_start, y_start, theta, phi, h = geometry.H, s_position = geometry.position)
    f = event_generator_functions.check_particle_position(x_stop, y_stop, l = geometry.L, w = geometry.W, s_position = geometry.position) 
    
    mask = numpy.full(N, True)
    
    data = numpy.vstack((x_start, y_start, theta, phi, x_stop, y_stop, f)).T    
    print("x_start, y_start, theta, phi, x_stop, x_stop, flag \n", data)                  
    epsilon = numpy.sum(f)/N    
    print("Number of events hitting the scintillator/Total number of events on S1:", numpy.sum(f), "/", N, "=", epsilon)

    if(output_file_events.endswith('.txt')): 
      header ='%s \nE[MeV], P [MeV], beta, x1[m], y1[m], theta, phi, x3[m], y3[m], flag\n' % datetime.datetime.now()
      fmt = ['%.4f', '%.4f', '%.4f', '%.4f', '%.6f', '%.2f', '%.2f', '%.4f', '%.6f', '%d']
      numpy.savetxt(output_file_events, numpy.transpose([E_kin[mask], P[mask], beta[mask], x_start[mask], y_start[mask], theta[mask], phi[mask], x_stop[mask], y_stop[mask], f[mask]]) , fmt=fmt, header=header)
      print("Output file saved!\n\n")   
    print("Time of execution: %s seconds " % (time.time() - start_time))
           

      
