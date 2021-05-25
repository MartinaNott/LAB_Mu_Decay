import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt


import plot_functions
import functions
import utilities


def efficiency(Nsel, Ntot): 
    epsilon = Nsel/Ntot
    epsilon_err = numpy.sqrt(epsilon * (1. - epsilon)/Ntot)
    epsilon = epsilon * 100
    epsilon_err = epsilon_err*100
    #print(epsilon, epsilon_err)
    return epsilon, epsilon_err
    

#Misure Martina    
N_sel = numpy.array([573, 573, 864, 581, 587, 587 ])
N_tot = numpy.array([1774, 603, 890, 617, 670, 1051 ])
epsilon_d1, epsilon_err_d1 = efficiency(N_sel, N_tot)    

#Misure 6 maggio
N_sel = numpy.array([2121, 2788, 2686, 2556, 2663, 2617 ])
N_tot = numpy.array([7247, 2893, 3674, 2640, 2898, 4351 ])
epsilon_d2, epsilon_err_d2 = efficiency(N_sel, N_tot)    

#Misure 11 maggio 
N_sel = numpy.array([856, 2293, 2127, 2063, 2129, 2063 ])
N_tot = numpy.array([2930, 2355, 2856, 2142, 2324, 3459 ])
epsilon_d3, epsilon_err_d3 = efficiency(N_sel, N_tot)    

#Misure 18 maggio 
N_sel = numpy.array([2502, 2140, 2904, 2083, 2114, 2450])
N_tot = numpy.array([9103, 2210, 3211, 2200, 2208, 3636])
epsilon_d4, epsilon_err_d4 = efficiency(N_sel, N_tot)    
    
days = numpy.array([1., 2., 3., 4.])


epsilon_matrix = numpy.vstack([epsilon_d1, epsilon_d2, epsilon_d3, epsilon_d4])
epsilon_err_matrix = numpy.vstack([epsilon_err_d1, epsilon_err_d2, epsilon_err_d3, epsilon_err_d4])

plt.figure()
for i in range(len(epsilon_d1)):
    legend = '%d' %(i+1)
    plot_functions.scatter_plot(days, epsilon_matrix[:, i]/100, 'days', 'epsilon', dx = None, dy = epsilon_err_matrix[:, i]/100,  title = '', fmt='o-', legend = legend)  
    epsilon_mean = numpy.mean(epsilon_matrix[:, i]/100)
    #epsilon_stdev = 
    print("epsilon_mean: ", epsilon_mean)
plt.ion()    
plt.show()    
