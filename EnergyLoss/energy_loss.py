import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

import constants 

LEAD_DENSITY = 11.34 # g/cm^3
ALUMINIUM_DENSITY = 2.7 # g/cm^3
IRON_DENSITY = 7.874 # g/cm^3
CARBON_DENSITY = 2. # g/cm^3


def inverse_stp_pwr(x):
    return 1./stp_pwr_spline(x)

def collision_range(initial_kin_energy):
    y, abserr = quad(inverse_stp_pwr, 0.01, initial_kin_energy)
    return y #g/cm^2
    
def collision_range_muon(initial_kin_energy):
    y, abserr = quad(inverse_stp_pwr, 2.1, initial_kin_energy )
    return y #g/cm^2

density_data_dict = {'lead' : ('estar_data_lead.txt', LEAD_DENSITY)}
                     #'aluminium' : ('estar_data_aluminium.txt', ALUMINIUM_DENSITY) ,
                     #'iron' : ('estar_data_iron.txt', IRON_DENSITY),
                     #'carbon' : ('estar_data_carbon.txt', CARBON_DENSITY)}

for element in density_data_dict:
    data_file = density_data_dict[element][0]
    density = density_data_dict[element][1]    
    kin_ene, collision_stp_pwr, rad_stop_pwr, total_stp_pwr, csda_range = numpy.loadtxt(data_file, unpack=True)
    csda_range = csda_range/density
    
    #PERDITA DI ENERGIA
    plt.figure()
    plt.title(element, fontsize=14)
    plt.plot(kin_ene, collision_stp_pwr, 'o')
    stp_pwr_spline = interp1d(kin_ene, collision_stp_pwr, kind='cubic')
    ene_grid = numpy.logspace(-2, 3., 201)
    plt.plot(ene_grid, stp_pwr_spline(ene_grid), '-b', label = 'electron')

    muon_kin_ene = kin_ene *  constants.MUON_MASS / constants.ELECTRON_MASS
    plt.plot(muon_kin_ene, collision_stp_pwr, '*')
    stp_pwr_spline_muon = interp1d(muon_kin_ene, collision_stp_pwr, kind='cubic')
    ene_grid_muon = numpy.logspace(numpy.log10(2.1), numpy.log10(20000), 201)
    plt.plot(ene_grid_muon, stp_pwr_spline_muon(ene_grid_muon), '-r', label = 'muon')

    plt.legend()
    plt.xscale('log')
    plt.xlabel('Kinetic energy [MeV]', fontsize=14)
    plt.ylabel('Collision stopping power [MeV cm2/g]', fontsize=14)
    plt.yticks(fontsize=14, rotation=0)
    plt.xticks(fontsize=14, rotation=0) 
    plt.subplots_adjust(bottom = 0.13, left = 0.15)  


    #RANGE
    stp_pwr_spline = interp1d(kin_ene, collision_stp_pwr, kind='cubic')
    est_range = numpy.full(len(kin_ene), 0.)
    for i in range(len(kin_ene)):
        est_range[i] = collision_range(kin_ene[i])
    est_range = est_range/density #cm


    stp_pwr_spline = interp1d(muon_kin_ene, collision_stp_pwr, kind='cubic')
    est_muon_range = numpy.full(len(muon_kin_ene), 0.)
    for i in range(len(muon_kin_ene)):
        est_muon_range[i] = collision_range_muon(muon_kin_ene[i] )
    est_muon_range = est_muon_range/density #cm

    plt.figure()
    plt.title(element, fontsize=14 )
    plt.plot(kin_ene, csda_range, 'r-', label = 'range e (total)')
    plt.plot(kin_ene, est_range, 'r--', label = 'range e (collision)')
    plt.plot(muon_kin_ene, est_muon_range, 'b-', label = 'range mu (collision)')
    plt.xscale('log')
    plt.xlabel("Kinetic energy [MeV]", fontsize=14)
    plt.ylabel("Range [cm]", fontsize=14)
    plt.yticks(fontsize=14, rotation=0)
    plt.xticks(fontsize=14, rotation=0) 
    plt.subplots_adjust(bottom = 0.13, left = 0.15)  
    plt.ylim(0., 20.)
    plt.xlim(1., 1000.)
    plt.legend()


    output_file = 'range_%s.txt' % element
    header ='kin_ene, range_e, muon_kin_ene, range_mu\n' 
    fmt = ['%.4f', '%.4f', '%.4f', '%.4f']
    numpy.savetxt(output_file, numpy.transpose([kin_ene, csda_range, muon_kin_ene, est_muon_range]) , fmt=fmt, header=header)
    print("Output file saved!\n\n")   



plt.ion()
plt.show()

