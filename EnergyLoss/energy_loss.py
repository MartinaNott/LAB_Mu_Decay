import numpy
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

def inverse_stp_pwr(x):
    return 1./stp_pwr_spline(x)

def collision_range(initial_kin_energy):
    y, abserr = quad(inverse_stp_pwr, 0.01, initial_kin_energy)
    return y #g/cm^2
    
def collision_range_muon(initial_kin_energy):
    y, abserr = quad(inverse_stp_pwr, 2.1, initial_kin_energy )
    return y #g/cm^2


data_file = 'estar_data_lead.txt'
kin_ene, collision_stp_pwr, rad_stop_pwr, total_stp_pwr, csda_range = numpy.loadtxt(data_file, unpack=True)

LEAD_DENSITY = 11.34 # g/cm^3
ALUMINIUM_DENSITY = 2.7 # g/cm^3
IRON_DENSITY = 7.874 # g/cm^3

MUON_MASS = 105.658 # MeV/c^2
ELECTRON_MASS = 0.511 # MeV/c^2


#ELETTRONI
plt.figure("electron_energy_loss")
plt.plot(kin_ene, collision_stp_pwr, 'o')
stp_pwr_spline = interp1d(kin_ene, collision_stp_pwr, kind='cubic')
ene_grid = numpy.logspace(-2, 3., 201)
print(ene_grid)
plt.plot(ene_grid, stp_pwr_spline(ene_grid))
plt.xscale('log')
plt.xlabel('kinetic energy [MeV]')
plt.ylabel('collision stopping power [MeV cm2/g]')

est_range = numpy.full(len(kin_ene), 0.)
for i in range(len(kin_ene)):
    est_range[i] = collision_range(kin_ene[i])
est_range = est_range/LEAD_DENSITY #cm

plt.figure("Electron_range")
plt.plot(kin_ene, csda_range/LEAD_DENSITY, 'r-')
plt.plot(kin_ene, est_range, 'b-')
plt.xscale('log')
plt.xlabel("kinetic energy [MeV]")
plt.ylabel("Range [cm]")
plt.ylim(0., 6.)


#MUONI

muon_kin_ene = kin_ene *  MUON_MASS / ELECTRON_MASS
plt.figure("muon_energy_loss")
plt.plot(muon_kin_ene, collision_stp_pwr, 'o')
stp_pwr_spline = interp1d(muon_kin_ene, collision_stp_pwr, kind='cubic')
ene_grid = numpy.logspace(numpy.log10(2.1), numpy.log10(20000), 201)
print(muon_kin_ene)
plt.plot(ene_grid, stp_pwr_spline(ene_grid))
plt.xscale('log')
plt.xlabel('kinetic energy [MeV]')
plt.ylabel('collision stopping power [MeV cm2/g]')

est_muon_range = numpy.full(len(muon_kin_ene), 0.)
for i in range(len(muon_kin_ene)):
    est_muon_range[i] = collision_range_muon(muon_kin_ene[i] )
est_muon_range = est_muon_range/LEAD_DENSITY #cm


plt.figure("Muon_range")
plt.plot(muon_kin_ene, est_muon_range, 'b-')
plt.xscale('log')
plt.xlabel("Energia cinetica [MeV]", fontsize=14)
plt.ylabel("Range [cm]", fontsize=14)
plt.ylim(0., 2.)
plt.xlim(0., 100.)

plt.yticks(fontsize=14, rotation=0)
plt.xticks(fontsize=14, rotation=0) 
plt.subplots_adjust(bottom = 0.13, left = 0.15)  

plt.ion()
plt.show()

