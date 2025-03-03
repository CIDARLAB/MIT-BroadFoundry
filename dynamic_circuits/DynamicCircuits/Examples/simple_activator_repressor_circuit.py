# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:30:30 2015

@author: Alex Lim
"""
#Circuit Models: Simple Activator and Repressor

# cell transcription circuit modeling
import numpy as np
import matplotlib.pyplot as plt
import inputs
from scipy.integrate import odeint
plt.ion()

#Simple cases where K = 1, n = 1 for Hill Equations
xMax = 200
yMax = 55
B = 1                       # rate of production of Y
a = 0.02                    # rate of degradation/dilution
Yst = B/a                   # steady state
Tr = np.log(2)/a            # response time
init = [0,Yst]              # initial concentrations
        
t = np.linspace(0, xMax, 1000)  # time grid (x-axis)


# solve the system of DEs. X is the input signal given by 'inputs' file
def f(ini, t):
        Xi = inputs.linInput(t,0,1)
        Yi = ini[0]
        Zi = ini[1]
        # the model equations
        f0 = B*Xi/(1+Xi) - a*Yi
        f1 = B/(1+Xi) - a*Zi
        return [f0,f1]

# solve the DE
soln = odeint(f, init, t)
X =inputs.linInput(t,0,10)
Y = soln[:,0]
Z = soln[:,1]

# Simple Activator, K=1, n=1
plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,Y)
plt.axhline(y=Yst, color='r', ls='dashed')
plt.title('Simple Activator')
plt.xlabel('Time')
plt.ylabel('Concentration Y')

# Simple Repressor, K=1, n=1
plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,Z)
plt.axhline(y=Yst, color='r', ls='dashed')
plt.title('Simple Repressor')
plt.xlabel('Time')
plt.ylabel('Concentration Z')
