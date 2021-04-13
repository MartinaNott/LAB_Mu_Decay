import sys
sys.path.insert(1, '/home/testaovo/Scrivania/LABORATORIO/muon_decay/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt

import plot_functions

date = '13/04/21'

data_dict = {'Scintillator1' : ('Vth_counts_s1_1700V.txt', 'Vth_counts_s1_1700V.txt' ), 
                     'Scintillator2' : ('Vth_counts_s2_1700V.txt', 'Vth_counts_s1_1700V.txt') ,
                     'Scintillator3' : ('Vth_counts_s3_1700V.txt', 'Vth_counts_s1_1700V.txt'),
                     'Scintillator4' : ('Vth_counts_s4_1700V.txt', 'Vth_counts_s1_1700V.txt'), 
                     'Scintillator5' : ('Vth_counts_s5_1700V.txt', 'Vth_counts_s1_1700V.txt'), 
                     'Scintillator6' : ('Vth_counts_s6_1700V.txt', 'Vth_counts_s1_1700V.txt') }

plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][0]
    Vth, counts, time = numpy.loadtxt(data_file, unpack=True)
    Vss = 1700. #V
    title = '%s - %d V e %f s' % (date, Vss, time[0])
             
    dcounts = numpy.sqrt(counts)
    plt.errorbar(Vth, counts, yerr = dcounts, label = scintillator)
    plot_functions.set_plot('$V_{th}$', 'Counts', title = title)

"""
"plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][1]
    Vss = #V
    title = '%s - %d V e %f s' % (date, Vss, time[0])
    Vth, counts, time = numpy.loadtxt(data_file, unpack=True)         
    dcounts = numpy.sqrt(counts)
    plt.errorbar(Vth, counts, yerr = dcounts, label = scintillator)
    plot_functions.set_plot('$V_{th}$', 'Counts', title = title)
"""


plt.ion()
plt.show()

