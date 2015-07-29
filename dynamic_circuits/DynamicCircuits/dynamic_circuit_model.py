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
n = 200.              # hill coefficient
Ky = 0.5            # activation coefficient of Y
Kz = 0.5            # repression coefficient of Z

t = np.linspace(xMin, xMax, itr)   # time grid

per = 4
amp = 1
dis = 1
X = inputs.sinInput(t,per,amp,0)
X2 = inputs.squInput(t,per,amp,dis)

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

plt.figure()
plt.axis([xMin,xMax,yMin,yMax])
plt.plot(t,X,'b',label = 'input')
plt.plot(t,X2,'b--',label = 'digital')
plt.plot(t,Y,'g',label = 'activator')
plt.plot(t,Z,'r',label = 'repressor')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input & Output')
plt.legend(loc=0)

'''
# Input Plot
plt.figure()
plt.axis([xMin,xMax,yMin,yMax])
plt.plot(t,X,'b')
plt.plot(t,X2,'b--')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input: X concentration')

# Plot of outputs: activity (activator and repressor)
plt.figure()
plt.axis([xMin,xMax,0,3])
plt.plot(t,Y, 'g', label = 'activator')
plt.plot(t,Z, 'r', label = 'repressor')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)
'''