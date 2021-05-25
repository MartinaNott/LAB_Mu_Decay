import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt

import plot_functions

data_dict = {'Scintillator1' : ('dati_15_04/Vss_40mv_1.txt', 'dati_15_04/Vss_70mv_1.txt', 'dati_20_04/s1.txt', 'dati_20_04/s1.txt' ), 
             'Scintillator2' : ('dati_15_04/Vss_40mv_2.txt', 'dati_15_04/Vss_70mv_2.txt', 'dati_20_04/s2.txt', 'dati_20_04/s2.txt') ,
             'Scintillator3' : ('dati_15_04/Vss_40mv_3.txt', 'dati_15_04/Vss_70mv_3.txt', 'dati_20_04/s3.txt', 'dati_20_04/s3_13.txt'),
             'Scintillator4' : ('dati_15_04/Vss_40mv_4.txt', 'dati_15_04/Vss_70mv_4.txt', 'dati_20_04/s4.txt', 'dati_20_04/s4_24.txt'), 
             'Scintillator5' : ('dati_15_04/Vss_40mv_5.txt', 'dati_15_04/Vss_70mv_5.txt', 'dati_20_04/s5.txt', 'dati_20_04/s5.txt'),
             'Scintillator6' : ('dati_15_04/Vss_40mv_6.txt', 'dati_15_04/Vss_70mv_6.txt','dati_20_04/s6.txt', 'dati_20_04/s6.txt') }

plt.figure()
date = '15/04/21'
for scintillator in data_dict:
    data_file = data_dict[scintillator][0]
    Vss, counts, time = numpy.loadtxt(data_file, unpack=True)
    #Vth = 40. #mV
    title = '%s' % (date)
    rate = counts/time
    dcounts = numpy.sqrt(counts)
    drate = dcounts/time
    plt.errorbar(Vss, rate, yerr = drate, label = scintillator)
    plot_functions.set_plot('$V_a$ [V]', 'Rate [Hz]', title = title)
   
plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][1]
    Vss, counts, time = numpy.loadtxt(data_file, unpack=True)
    #Vth = 70. #mV
    title = '%s' % (date)
    rate = counts/time
    dcounts = numpy.sqrt(counts)
    drate = dcounts/time
    plt.errorbar(Vss, rate, yerr = drate, label = scintillator)
    plot_functions.set_plot('$V_a$ [V]', 'Rate [Hz]', title = title)

date = '20/04/21'
plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][2]
    Vss, counts, time = numpy.loadtxt(data_file, unpack=True)
    #Vth = 40. #mV
    title = '%s' % (date)
    rate = counts/time
    dcounts = numpy.sqrt(counts)
    drate = dcounts/time
    plt.errorbar(Vss, rate, yerr = drate, label = scintillator)
    plot_functions.set_plot('$V_a$ [V]', 'Rate [Hz]', title = title)

plt.figure()
for scintillator in data_dict:
    data_file = data_dict[scintillator][3]
    Vss, counts, time = numpy.loadtxt(data_file, unpack=True)
    #Vth = 40. #mV
    title = '%s' % (date)
    rate = counts/time
    dcounts = numpy.sqrt(counts)
    drate = dcounts/time
    plt.errorbar(Vss, rate, yerr = drate, label = scintillator)
    plot_functions.set_plot('$V_a$ [V]', 'Rate [Hz]', title = title)
            
plt.ion()
plt.show()

