# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 2:26:32 2015

@author: Alex Lim
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

xMax = 8
yMax = 1.1 
per = 4     # period of sinusoid
itr = 1000  # time iterations
y = [0]     # initial conditions
a = 0.2       # degradation rate
By = 1     # rate of production of Y
n = 1000      # hill coefficient
Ky = 0.5      # activation coefficient

        
t = np.linspace(0, xMax, itr)   # time grid
X = 0.5-0.5*np.cos(np.pi*t/(0.5*per))
B = By * np.power(X,n)/(np.power(X,n)+np.power(Ky,n))

y = [0]

# solve the system dx/dt = f(x, t)
def f(y, t):
        Yi = y[0]
        # the model equations
        f0 = By - a*Yi
        return [f0]

# solve the DEs
ySoln = odeint(f, y, t)
Y = ySoln[:,0]


plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,X)
plt.plot(t,B)

plt.figure()
plt.axis([0, 8, 0, 5])
plt.plot(t,Y)