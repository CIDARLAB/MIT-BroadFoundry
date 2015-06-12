# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:08:38 2015

@author: Alex Lim
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
plt.ion()

#Autoregulator
xMax = 200.
yMax = 100
B = 1.                       # rate of production of Y
B0 = 0.021                   # rate of auto-activation
a = 0.02                     # rate of degradation/dilution
X0 = [0,0,0]                 # initial concentration of Y
Yst=B/a                      # steady state
K = 50.                      # repression threshold
Tr = math.log(2)/a           # response time
n = 4                        # hill coefficient

# solve the system dy/dt = f(y, t)
def f(x, t):
        Xi = x[0]
        Yi = x[1]
        Zi = x[2]
        # the model equations
        f0 = B/(1+math.pow(Xi/K,n)) - a*Xi
        f1 = B - a*Yi
        f2 = B/(1+math.pow(Zi/K,n)) + (B0-a)*Zi
        return [f0,f1,f2]
        
t = np.linspace(0, xMax, 1000)   # time grid

# solve the DEs
soln = odeint(f, X0, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]

plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,Z, label = 'pos. autoregulator')
plt.plot(t,Y, label = 'simple regulator')
plt.plot(t,X, label = 'neg. autoregulator')
plt.axhline(y=K/2, xmax=Tr/xMax, color='r', ls='dashed')
plt.axvline(x=Tr, ymax=K/(2*yMax), color='r', ls='dashed')
plt.title('Circuit Comparisons')
plt.xlabel('Time')
plt.ylabel('[Concentration]')
plt.legend(loc=0)

plt.figure()
plt.axis([0,2,0,1])
plt.axhline(y=0.5, xmax=0.7, color='r', ls='dashed')
plt.axvline(x=0.675, ymax=0.5, color = 'r', ls='dashed')
plt.axvline(x=1, ymax=0.5, color = 'r', ls='dashed')
plt.axvline(x=1.4, ymax=0.5, color = 'r', ls='dashed')
plt.plot(t/Tr,Z/Z[999], label = 'pos. autoregulator')
plt.plot(t/Tr,Y/Yst, label = 'simple regulator')
plt.plot(t/Tr,X/X[999], label = 'neg. autoregulator')
plt.title('Circuit Comparisons')
plt.xlabel('t/Tr')
plt.ylabel('[X]/Xst')
plt.legend(loc=0)


