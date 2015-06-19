# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:30:30 2015

@author: Alex Lim
"""
#Circuit Models #1: Simple Activator and Repressor

# cell transcription circuit modeling
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
plt.ion()

# solve the system dy/dt = f(y, t)
def f(y, t):
        Yi = y[0]
        # the model equations
        f0 = B - a*Yi
        return f0

#Simple Activator
xMax = 200
yMax = 120
B = 1                       # rate of production of Y
a = 0.02                    # rate of degradation/dilution
Yst = B/a                   # steady state
Tr = math.log(2)/a          # response time
Y0 = [0]                    # initial concentration of Y
        
t = np.linspace(0, xMax, 1000)  # time grid (x-axis)
c = np.linspace(0, yMax, 1000)  # concentration grid (y-axis)
    
# solve the DE
soln = odeint(f, Y0, t)
Y = soln[:,0]

plt.figure()
plt.plot(t/Tr,Y/Yst)
plt.axhline(y=0.5, xmin=0, xmax=1./6, color='r', ls='--')
plt.axvline(x=1, ymin=0, ymax=0.5, color='r', ls='--')
plt.xlabel('t/Tr')
plt.ylabel('Y/Yst')
plt.title("Activation")

plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,Y, label = "Yst*(1-e^(-a*t))")
plt.plot(c, B*c, label = "B*t")
plt.axhline(y=Yst/2, xmin=0, xmax=Tr/xMax, color='r', ls='--')
plt.axhline(y=Yst, color='r', ls='dashed')
plt.axvline(x=Tr, ymin=0, ymax=Yst/(2*yMax), color='r', ls='--')
plt.title('Simple Activator')
plt.xlabel('Time')
plt.ylabel('Concentration Y')
plt.legend(loc=0)


# Simple Repressor, degradation, reduced production
B = 0                       # rate of production of Y, drops to 0
a = 0.02                    # rate of degradation/dilution
Tr = math.log(2)/a          # response time
Y0 = [Yst]                  # initial concentration of Y
        
# solve the DE
soln = odeint(f, Y0, t)
Y = soln[:]

plt.figure()
plt.plot(t/Tr,Y/Yst)
plt.axhline(y=0.5, xmin=0, xmax=1./6, color='r', ls='--')
plt.axvline(x=1, ymin=0, ymax=0.5, color='r', ls='--')
plt.xlabel('t/Tr')
plt.ylabel('Y/Yst')
plt.title("Repression")

plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,Y, label = "Yst*e^(-a*t)")
plt.axhline(y=Yst/2, xmin=0, xmax=Tr/xMax, color='r', ls='--')
plt.axvline(x=Tr, ymin=0, ymax=Yst/(2*yMax), color='r', ls='--')
plt.title('Simple Repressor')
plt.xlabel('Time')
plt.ylabel('Concentration Y')
plt.legend(loc=0)