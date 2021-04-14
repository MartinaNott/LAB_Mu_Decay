import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import geometry


def angle_generator(num_events, pdf, delta_2pi = False):      
  if delta_2pi is True: 
    theta = numpy.linspace(-numpy.pi, numpy.pi, 200)    
  else:
    theta = numpy.linspace(-numpy.pi/2, numpy.pi/2, 200) 
  cdf_y  = []
  for i in range(len(theta)):
    y, rest = quad(pdf, theta[0], theta[i])  
    cdf_y.append(y)
  cdf_y, unique_indices = numpy.unique(cdf_y, return_index=True)
  cdf_y /= cdf_y[-1]
  theta = theta[unique_indices] 
  funzione = interp1d(cdf_y, theta)       
  x = numpy.random.uniform(0., 1., num_events)
  theta_particle = funzione(x)
  phi_particle = numpy.random.uniform(0, 2 * numpy.pi, num_events )   
  return theta_particle, phi_particle

def energy_generator(num_events, pdf, e_min, e_max, *args): 
  e = numpy.linspace(e_min, e_max, 1000)  
  cdf_y  = numpy.full(len(e), 0.)
  for i in range(len(e)):
    y, rest = quad(pdf, e[0], e[i], args = args)
    cdf_y[i] = y       
  cdf_y, unique_indices = numpy.unique(cdf_y, return_index=True)
  cdf_y = cdf_y/cdf_y[-1]
  e = e[unique_indices]  
  ppf_spline = interp1d(cdf_y, e)         
  x = numpy.random.uniform(0., 1., num_events)
  E_kin = ppf_spline(x)  
  return E_kin

def p_relation(E_kin , mass ): 
  E_particle = E_kin + mass
  P_particle = numpy.sqrt( E_particle**2 - mass **2)
  return P_particle

def beta(P_particle, E_particle): 
  beta_particle = P_particle / E_particle  
  return beta_particle
  
def gamma(beta = -1, E_particle = -1 , mass = -1):   
  if (beta != -1) :
    gamma = numpy.sqrt(1/(1 - beta**2))  
    return gamma 
  if (E_particle != -1):
    gamma = E_particle/mass 
    return gamma   
    
def position_on_scintillator_generator(num_events, l = geometry.L, w = geometry.W): 
  x_s = numpy.random.uniform(0., l, num_events)
  y_s = numpy.random.uniform(0., w, num_events)
  return x_s, y_s
  
def cosmic_rays_propagation(x_s, y_s, theta, phi, h = geometry.H, s_position = geometry.position):
  z = h + s_position[2]
  x_new_s = x_s + numpy.cos(phi) * numpy.tan(theta) * z 
  y_new_s = y_s + numpy.sin(phi) * numpy.tan(theta) * z   
  return x_new_s, y_new_s, z
  
  
def check_particle_position(x_s, y_s, l = geometry.L, w = geometry.W, s_position = geometry.position):   
  mask_x = (x_s > s_position[0]) * (x_s < (s_position[0] + l))
  mask_y =  (y_s > s_position[1]) * (y_s < (s_position[1] + w)) 
  mask = mask_x * mask_y             
  return mask

def energy_range_functions(data_file): 
  electron_kin_ene, electron_ranges_in_material, muon_kin_ene, muon_ranges_in_material = numpy.loadtxt(data_file, unpack=True)  
  range_spline_electron = interp1d(electron_kin_ene, electron_ranges_in_material, kind='cubic')
  range_spline_muon = interp1d( muon_kin_ene, muon_ranges_in_material, kind='cubic')
  return range_spline_muon, range_spline_electron
  
def propagation_in_material(range_lenght, theta, phi, x_start = 0., y_start = 0., z_start = 0.): 
  x_stop = x_start + range_lenght * numpy.sin(theta) * numpy.cos(phi) 
  y_stop = y_start + range_lenght * numpy.sin(theta) * numpy.sin(phi)
  z_stop = z_start + range_lenght * numpy.cos(theta)  
  return x_stop, y_stop, z_stop
  
