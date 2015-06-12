# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:36:54 2015

@author: Alex Lim
"""

# cell transcription circuit modeling
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
plt.ion()

#Autoregulator
xMax = 160.
yMax = 60.
B = 1.                       # rate of production of Y
a = 0.02                    # rate of degradation/dilution
X0 = [0]                    # initial concentration of Y
K = 50.                     # repression threshold
Tr = K/(2*B)                  # response time
n = 1000000.                 # hill coefficient

# solve the system dy/dt = f(y, t)
def f(x, t):
        Xi = x[0]
        # the model equations
        f0 = B
        return f0
        
t = np.linspace(0, 2*Tr, 1000)   # time grid
t2 = np.linspace(0, 160, 1000)
# solve the DEs
soln = odeint(f, X0, t)
soln2 = odeint(f, X0, t2)
X = soln[:,0]
X2 = soln2[:,0]

plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,X)
plt.plot(t2,X2, ls='dashed')
plt.axhline(y=K, xmin=2*Tr/xMax)
plt.axhline(y=K/2, xmax=Tr/xMax, color='r', ls='dashed')
plt.axvline(x=Tr, ymax=K/(2*yMax), color='r', ls='dashed')
plt.title('Autoregulator')
plt.xlabel('Time')
plt.ylabel('Concentration X')
plt.legend(loc=0)





#AND-Gate C1-FFL

#OR-Gate C1-FFL
