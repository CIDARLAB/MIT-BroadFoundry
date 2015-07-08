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
yMax = 3
itr = 1000           # time iterations
inic = [0,0]        # initial conditions
By = 1.             # rate of production of Y
ay = 0.5            # degradation rate of Y
Bz = 1.             # rate of production of Z
az = 0.5            # degradation rate of Z
n = 2.              # hill coefficient
Ky = 0.5            # activation coefficient of Y
Kz = 0.5            # activation coefficient of Z

t = np.linspace(xMin, xMax, itr)   # time grid

X = inputs.sinInput(t)

# solve the system dx/dt = f(x, t)
def f(ic, t):
        Xi = inputs.sinInput(t)
        Yi = ic[0]
        Zi = ic[1]    
        # Activator Equation
        f0 = By*(np.power(Xi,n)/(np.power(Xi,n)+np.power(Ky,n))) - ay*Yi
        # Repressor Equation        
        f1 = Bz/(1+np.power(Xi/Kz,n)) - az*Zi
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
plt.axis([xMin,xMax,0,2])
plt.plot(t,Y, label = 'activator')
plt.plot(t,Z, label = 'repressor')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)

plt.figure()
plt.axis([xMin,xMax,yMin,yMax])
plt.plot(t,X, label = 'input')
plt.plot(t,Y, label = 'X=activator')
plt.plot(t,Z, label = 'X=repressor')
plt.axvline(x=1, ls='--')
plt.axvline(x=2, ls='--')
plt.axvline(x=3, ls='--')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input & Output')
plt.legend(loc=0)