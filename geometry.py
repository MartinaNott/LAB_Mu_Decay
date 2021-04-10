import numpy

c = 29979245800  #cm/s

#Large scintillators dimensions
L = 50. #cm
W = 40. #cm
H = 1.  #cm

#Distance between scintillators: 
h_ij = 5. #cm


#Little scintillator dimensions: 
L_little = 30.  #cm
W_little = 30.  #cm
Z_little = 10.  #cm


position = numpy.array([0., 0., h_ij])
