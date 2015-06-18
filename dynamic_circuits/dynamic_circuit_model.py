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

xMin = 0
xMax = 10
yMin = 0
yMax = 4
itr = 100          # time iterations
inic = [0,0]        # initial conditions
By = 1.             # rate of production of Y
ay = 0.5            # degradation rate of Y
Bz = 1.             # rate of production of Z
az = 0.5            # degradation rate of Z
n = 1.              # hill coefficient
Ky = 1.             # activation coefficient
Ky = 0.2             # activation coefficient

t = np.linspace(0, xMax, itr)   # time grid
    

X = inputs.squInput(t)

# solve the system dx/dt = f(x, t)
def f(y, t):
        Xi = inputs.squInput(t)
        Yi = y[0]
        Zi = y[1]    
        # Activator Equation
        f0 = By*(np.power(Xi,n)/(np.power(Xi,n)+np.power(Ky,n))) - ay*Yi
        # Repressor Equation        
        f1 = Bz/(1+np.power(Xi/Ky,n)) - az*Zi
        return [f0 , f1]

# solve the DEs
soln = odeint(f, inic, t)
Y = soln[:,0]
Z = soln[:,1]

# Input Plot
plt.figure()
plt.axis([xMin,xMax,yMin,yMax])
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

