import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import matplotlib

import numpy
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import event_generator_functions
import geometry
import functions
import constants
import plot_functions


N = 1000000
z_material = (0., 30.)
n_bin = 300
material_depth = numpy.linspace(*z_material, n_bin) 

element_dict = {'lead' : '../EnergyLoss/range_lead.txt',
                     'aluminium' : '../EnergyLoss/range_aluminium.txt' ,
                     'iron' : '../EnergyLoss/range_iron.txt',
                     'carbon' : '../EnergyLoss/range_carbon.txt'}
                     
#E_kin_mu = event_generator_functions.energy_generator(N, functions.spectrum, 1., 1.e5, -2.7, 300., 50.) 
E_kin_mu = event_generator_functions.energy_generator(N, functions.muon_spectrum, 1., 1000) 
#E_kin_mu = numpy.random.uniform(3, 100, N )   
#plt.figure()
#plt.hist(E_kin_mu, range = (0., 1000), bins = 100)

theta_mu, phi_mu = event_generator_functions.angle_generator(N, functions.dist_theta)   


for element in element_dict:
    data_file = element_dict[element]
    spline_muon, spline_electron = event_generator_functions.energy_range_functions(data_file)
    range_mu = spline_muon(E_kin_mu) 
    epsilon_muon = []
    epsilon_electron = []

    for i in range (0, n_bin):       
        z_range_mu = range_mu * numpy.cos(theta_mu)
        mu_stopped_mask = (z_range_mu < material_depth[i])
        mu_stopped = numpy.sum(mu_stopped_mask)  
        epsilon_muon.append(mu_stopped/N)
        z_range_mu = z_range_mu[mu_stopped_mask]  
            
        if mu_stopped != 0: 
            E_kin_ele = event_generator_functions.energy_generator(mu_stopped, functions.electron_spectrum, 0., constants.MUON_MASS * 0.5)
            range_ele = spline_electron(E_kin_ele)  
            theta_ele, phi_ele = event_generator_functions.angle_generator(mu_stopped, functions.dist_theta, delta_2pi = True)   
            z_abs_ele = z_range_mu + range_ele * numpy.cos(theta_ele)
            ele_stopped_mask = (0. < z_abs_ele) * (z_abs_ele < material_depth[i])
            ele_stopped = numpy.sum(ele_stopped_mask)  
            epsilon_electron.append(1. - ele_stopped/mu_stopped)  
        else:
            epsilon_electron.append(0.)

    material_depth = numpy.array(material_depth)   
    epsilon_electron = numpy.array(epsilon_electron)   
    epsilon_muon = numpy.array(epsilon_muon) 
    epsilon = epsilon_electron * epsilon_muon
    
    
    label_ele = 'e in %s' % element
    label_mu = 'mu in %s' % element
    plt.figure()
    plt.plot(material_depth, epsilon_electron, '-r', label = label_ele )
    plt.plot(material_depth, epsilon_muon, '-b', label = label_mu)
    plot_functions.set_plot("material_depth [cm]", "N_stopped/N_tot", title = element )
    
    plt.figure()
    plt.plot(material_depth, epsilon, '-')
    plot_functions.set_plot("material_depth [cm]", "epsilon", title = element )
    

plt.ion()
plt.show()


