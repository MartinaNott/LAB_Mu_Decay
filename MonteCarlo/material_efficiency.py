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

save_fig = False
N = 1000000
z_material = (0., 80.)
n_bin = 100
material_depth = numpy.linspace(*z_material, n_bin) 

element_dict = {'lead' : '../EnergyLoss/range_lead.txt',
                'aluminium' : '../EnergyLoss/range_aluminium.txt' ,
                'iron' : '../EnergyLoss/range_iron.txt',
                'carbon' : '../EnergyLoss/range_carbon.txt'}

E_min = 10
E_max = 1000        
#E_kin_mu = event_generator_functions.energy_generator(N, functions.spectrum, 1., 1.e5, -2.7, 300., 50.) 
E_kin_mu = numpy.random.uniform(E_min, E_max, N )   

#plt.figure()
#plot_functions.plot_histogram(E_kin_mu, "$", "", bins = 100, range = , title = '')

theta_mu, phi_mu = event_generator_functions.angle_generator(N, functions.dist_theta)   

for element in element_dict:
    data_file = element_dict[element]
    spline_muon, spline_electron = event_generator_functions.energy_range_functions(data_file)
    range_mu = spline_muon(E_kin_mu) 
    epsilon_muon = []
    epsilon_electron = []
    electron_up = []
    electron_down = []

    for i in range (0, n_bin):       
        z_range_mu = range_mu * numpy.cos(theta_mu)
        mu_stopped_mask = (z_range_mu > 0.) * (z_range_mu < material_depth[i].max())
        mu_stopped = numpy.sum(mu_stopped_mask)  
        epsilon_muon.append(mu_stopped/N)
        z_range_mu = z_range_mu[mu_stopped_mask]  
        
        #plt.figure()
        #title = 'muon %f' % material_depth[i]
        #plot_functions.plot_histogram(z_range_mu, "$z_{muon}$ [cm]", "", bins = 100, range = z_material, title = title)
             
        if mu_stopped != 0: 
            E_kin_ele = event_generator_functions.energy_generator(mu_stopped, functions.electron_spectrum, 0., constants.MUON_MASS * 0.5)
            range_ele = spline_electron(E_kin_ele)  
            theta_ele, phi_ele = event_generator_functions.angle_generator(mu_stopped, functions.dist_theta, delta_2pi = True)   
            z_abs_ele = z_range_mu + range_ele * numpy.cos(theta_ele)
            ele_passed_mask_up = (0. > z_abs_ele)
            ele_passed_mask_down = (z_abs_ele > material_depth[i])
            ele_stopped_mask = (0. < z_abs_ele) * (z_abs_ele < material_depth[i])
            ele_passed_up = numpy.sum(ele_passed_mask_up) 
            ele_passed_down = numpy.sum(ele_passed_mask_down) 
            electron_up.append(ele_passed_up/mu_stopped)  
            electron_down.append(ele_passed_down/mu_stopped)  
            ele_stopped = numpy.sum(ele_stopped_mask)  
            epsilon_electron.append(1. - ele_stopped/mu_stopped)  

            #plt.figure()
            #title = 'electron %f' % material_depth[i]
            #plot_functions.plot_histogram(z_abs_ele[~ele_stopped_mask], "$z_{ele}$ [cm]", "", bins = 100, title = title)
          
        else:
            epsilon_electron.append(0.)
            electron_up.append(-25.)
            electron_down.append(-25.)
            
    material_depth = numpy.array(material_depth)   
    electron_up = numpy.array(electron_up)
    electron_down = numpy.array(electron_down)
    epsilon_electron = numpy.array(epsilon_electron)   
    epsilon_muon = numpy.array(epsilon_muon) 
    epsilon = epsilon_electron * epsilon_muon
    

    title = '%s sp. piatto (%d-%d Mev)' % (element, E_min, E_max)    
    plt.figure()
    plt.plot(material_depth, epsilon_electron, '-r', label = 'electron' )
    plt.plot(material_depth, epsilon_muon, '-b', label = 'muon')
    plot_functions.set_plot("material_depth [cm]", "$N_i$/$N_{tot}$", title = title )
    
    if save_fig ==True:    
      plt.savefig('mu_ele_epsilon_%s.pdf' % element, format = 'pdf')    
    
    plt.figure()
    plt.plot(material_depth, epsilon, '-')
    plot_functions.set_plot("material_depth [cm]", "$\epsilon$", title = title )
    
    if save_fig ==True:    
      plt.savefig('epsilon_%s.pdf' % element, format = 'pdf')    
     
    plt.figure()
    plt.plot(material_depth, electron_up, '-', label = 'up')
    plt.plot(material_depth, electron_down, '-', label = 'down')
    plt.ylim(-0.1, 1.1)
    plot_functions.set_plot("material_depth [cm]", "$\$", title = title )
    
    #if save_fig ==True:    
    #  plt.savefig('epsilon_%s.pdf' % element, format = 'pdf')    



plt.ion()
plt.show()


