import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt

import plot_functions

date = '15/04/21'

data_dict = {'Scintillator1' : ('Vss_40mv_1.txt', 'Vss_1.txt' ), 
             'Scintillator2' : ('Vss_40mv_2.txt', 'Vss_2.txt') ,
             'Scintillator3' : ('Vss_40mv_3.txt', 'Vss_3.txt'),
             'Scintillator4' : ('Vss_40mv_4.txt', 'Vss_4.txt'), 
             'Scintillator5' : ('Vss_40mv_5.txt', 'Vss_5.txt'),
             'Scintillator6' : ('Vss_40mv_6.txt', 'Vss_6.txt') }

plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][0]
    Vss, counts, time = numpy.loadtxt(data_file, unpack=True)
    Vth = 40. #mV
    title = '%s - %d V ' % (date, Vth)
    rate = counts/time
    dcounts = numpy.sqrt(counts)
    drate = dcounts/time
    plt.errorbar(Vss, rate, yerr = drate, label = scintillator)
    plot_functions.set_plot('$V_{th}$ [mV]', 'Rate [Hz]', title = title)
   
plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][1]
    Vss, counts, time = numpy.loadtxt(data_file, unpack=True)
    Vth = 70. #mV
    title = '%s - %d V ' % (date, Vth)
    rate = counts/time
    dcounts = numpy.sqrt(counts)
    drate = dcounts/time
    plt.errorbar(Vss, rate, yerr = drate, label = scintillator)
    plot_functions.set_plot('$V_{th}$ [mV]', 'Rate [Hz]', title = title)
   



plt.ion()
plt.show()

