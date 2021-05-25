import sys
sys.path.insert(1, '/home/ele/lab/LAB_Mu_Decay')

import numpy
import matplotlib.pyplot as plt


import plot_functions
import functions
import utilities


N_1 = numpy.array([29592/291, 23514/60, 27370/377])
N_2 = numpy.array([383850/385, 299636/324, 221795/306])
N_3 = numpy.array([39937/152, 42673/120, 44482/166])
N_4 = numpy.array([701930/151, 615246/120, 527184/100, 62384/100, 98885/126])
N_5 = numpy.array([38041/151, 31289/122, 29685/130])
N_6 = numpy.array([18258/151, 14259/120, 18306/153]) 
dN_1 = numpy.sqrt(N_1)
dN_2 = numpy.sqrt(N_2)
dN_3 = numpy.sqrt(N_3)
dN_4 = numpy.sqrt(N_4)
dN_5 = numpy.sqrt(N_5)
dN_6 = numpy.sqrt(N_6)

days = numpy.array([1., 2., 3.])
days_4 = numpy.array([1., 2., 2.5, 2.6, 3.])

plt.figure()
plot_functions.scatter_plot(days, N_1, 'days', ' Noise[Hz]', dx = None, dy = dN_1,  title = '', fmt='o-', legend = 'S1')  
plot_functions.scatter_plot(days, N_2, 'days', ' Noise[Hz]', dx = None, dy = dN_2,  title = '', fmt='o-', legend = 'S2')  
plot_functions.scatter_plot(days, N_3, 'days', ' Noise[Hz]', dx = None, dy = dN_3,  title = '', fmt='o-', legend = 'S3')  
plot_functions.scatter_plot(days_4, N_4, 'days', ' Noise[Hz]', dx = None, dy = dN_4,  title = '', fmt='o-', legend = 'S4')  
plot_functions.scatter_plot(days, N_5, 'days', ' Noise[Hz]', dx = None, dy = dN_5,  title = '', fmt='o-', legend = 'S5')  
plot_functions.scatter_plot(days, N_6, 'days', ' Noise[Hz]', dx = None, dy = dN_6,  title = '', fmt='o-', legend = 'S6')  

plt.ion()    
plt.show()    
