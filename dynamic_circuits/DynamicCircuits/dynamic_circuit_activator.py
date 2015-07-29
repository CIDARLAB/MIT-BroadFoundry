# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:34 2015

@author: Alex Lim
"""

import matplotlib.pyplot as plt
import numpy as np
import inputs
from scipy.integrate import odeint
plt.ion()

xMax = 10
yMax = 4
itr = 500           # time iterations/resolution
inic = [1,0]        # initial conditions
By = 1.             # rate of production of Y
ay = 0.5            # degradation rate of Y
Bz = 1.             # rate of production of Z
az = 0.5            # degradation rate of Z
n = 2.              # hill coefficient
K = 1.              # activation coefficient
st = 2              # starting x-position
dis = 1             # distance
mag = 2             # magnitude

t = np.linspace(0, xMax, itr)   # time grid
X = inputs.stepInput(t, st, dis, mag)

# solve the system dx/dt = f(x, t)
def f(y, t):
        Xi = inputs.stepInput(t, st, dis, mag)
        Yi = y[0]
        Zi = y[1]    
        # Activator Equation
        f0 = By*(np.power(Xi,n)/(np.power(Xi,n)+np.power(K,n))) - ay*Yi
        # Repressor Equation        
        f1 = Bz/(1+np.power(Xi/K,n)) - az*Zi
        return [f0 , f1]

# solve the DEs
soln = odeint(f, inic, t)
Y = soln[:,0]
Z = soln[:,1]

# Input Plot: X-->Y and X--|Z
plt.figure()
plt.axis([0,xMax,0,yMax])
plt.plot(t,X)
plt.axhline(y=K, color='r', ls='--')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input: X concentration')

# Plot of outputs: Hill Function Response activity (activator and repressor)
plt.figure()
plt.axis([0,xMax,0,2])
plt.plot(t,Y, label = 'activator')
plt.plot(t,Z, color='r', label = 'repressor')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Hill Function Response')
plt.legend(loc=0)

