# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:34 2015

@author: Alex Lim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

per = 4     # period of sinusoid 
amp = 2     # amplitude of signal
nwaves = 2

xMax = per*nwaves
yMax = amp+0.1

itr = 1000  # time iterations
inic = [0,0,0]     # initial conditions
Bx = 1     # rate of production of X
ax = 0.2   # degradation rate of X
By = 1     # rate of production of Y
ay = 0.5      # degradation rate of Y
Bz = 1     # rate of production of Z
az = 0.5    # degradation rate of Z
n = 4     # hill coefficient
Ky = 1      # activation coefficient
       
t = np.linspace(0, xMax, itr)   # time grid

# solve the system dx/dt = f(x, t)
def f(inic, t):
        Xi = inic[0]
        Yi = inic[1]
        Zi = inic[2]
        # Input Signal Equation: wave function, ~.5*(amp-amp*cos(pi*t/(.5*per)))
        f0 = amp*np.pi/per*np.sin(np.pi*t/(0.5*per))
        # Activator Equation
        f1 = By*(np.power(Xi,n)/(np.power(Xi,n)+np.power(Ky,n))) - ay*Yi
        # Repressor Equation        
        f2 = Bz/(1+np.power(Xi/Ky,n)) - az*Zi
        return [f0 , f1, f2]

# solve the DEs
soln = odeint(f, inic, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]

# Input Plot
plt.figure()
plt.axis([0,xMax,-0.1,yMax])
plt.plot(t,X)
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input: X concentration')

# Plot of outputs: activity (activator and repressor)
plt.figure()
plt.axis([0,xMax,0,2])
plt.plot(t,Y, label = 'activator')
plt.plot(t,Z, label = 'repressor')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)

