# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:34 2015

@author: Alex Lim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

xMax = 10
yMax = 1.1 
per = 4     # period of sinusoid
itr = 1000  # time iterations
y = [0,0,0]     # initial conditions
ax = 0.
ay = 0.2       # degradation rate
Bx = 1     # rate of production of X
By = 1     # rate of production of Y
n = 1000     # hill coefficient
Ky = 5      # activation coefficient
       
t = np.linspace(0, xMax, itr)   # time grid

# solve the system dx/dt = f(x, t)
def f(y, t):
        Xi = y[0]
        Yi = y[1]
        Y2i = y[2]
        # the model equations
        f0 = np.pi/per*np.sin(np.pi*t/(0.5*per))
        f1 = By*(np.power(Xi,n)/(np.power(Xi,n)+np.power(Ky,n))) - ay*Yi
        f2 = By/(1+np.power(Xi/Ky,n)) - ay*Y2i
        return [f0 , f1, f2]

# solve the DEs
soln = odeint(f, y, t)
X = soln[:,0]
Y = soln[:,1]
Y2 = soln[:,2]

plt.figure()
plt.axis([0,xMax,-0.1,yMax])
plt.plot(t,X)
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input: X concentration')

plt.figure()
plt.plot(t,Y)
plt.plot(t,Y2)
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Y Activity')

