import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')
import numpy
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.interpolate import interp1d

import geometry

MUON_MASS = 105.658 #MeV

def cosmic_rays_angle_generator(num_events, pdf):      
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
  theta_muon = funzione(x)
  phi_muon = numpy.random.uniform(0, 2 * numpy.pi, num_events )   
  return theta_muon, phi_muon

def cosmic_rays_energy_generator(num_events, pdf, e_min, e_max, *args): 
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
  E_muon = E_kin + MUON_MASS
  P_muon = numpy.sqrt( E_muon**2 - MUON_MASS **2)
  beta_muon = P_muon / E_muon
  return E_muon, P_muon, beta_muon

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


