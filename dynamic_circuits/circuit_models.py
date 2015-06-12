# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:30:30 2015

@author: Alex Lim
"""


# cell transcription circuit modeling
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
plt.ion()

#Simple Regulator
#activation
B = 1                       # rate of production of Y
a = 0.02                    # rate of degradation/dilution
Yst = B/a                   # steady state
Tr = math.log(2)/a          # response time
Y0 = [0]                    # initial concentration of Y

# solve the system dy/dt = f(y, t)
def f(y, t):
        Yi = y[0]
        # the model equations
        f0 = B - a*Yi
        return f0
        
t = np.linspace(0, 160., 1000)   # time grid
t2 = np.linspace(0, 60., 1000)
    
# solve the DEs
soln = odeint(f, Y0, t)
Y = soln[:]

plt.figure()
plt.plot(t/Tr,Y/Yst)
plt.xlabel('t/Tr')
plt.ylabel('Y/Yst')
plt.legend(loc=0)

plt.figure()
plt.axis([0, 160, 0, 60])
plt.plot(t,Y, label = "Yst*(1-e^(-a*t))")
plt.plot(t2, B*t2, label = "B*t")
plt.axhline(y=Yst/2, xmin=0, xmax=Tr/160, color='r', ls='dashed')
plt.axhline(y=Yst, xmin=0, color='r', ls='dashed')
plt.axvline(x=Tr, ymin=0, ymax=Yst/(2*60), color='r', ls='dashed')
plt.title('Simple Regulator')
plt.xlabel('Time')
plt.ylabel('Concentration Y')
plt.legend(loc=0)


#degradation
B = 1                       # rate of production of Y
a = 0.02                    # rate of degradation/dilution
Yst = B/a                   # steady state
Tr = math.log(2)/a          # response time
Y0 = [0]                    # initial concentration of Y
def f(y, t):
        Yi = y[0]
        # the model equations
        f0 = B - a*Yi
        return f0

#Autoregulator
B = 0                       # rate of production of Y
a = 0.02                    # rate of degradation/dilution
Y0 = [0]                    # initial concentration of Y
K = 1                        # repression threshold
Yst = K                     # steady state
Tr = math.log(2)/a          # response time


#AND-Gate C1-FFL

#OR-Gate C1-FFL