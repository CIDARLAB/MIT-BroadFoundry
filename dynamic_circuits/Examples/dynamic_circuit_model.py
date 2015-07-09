# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:34 2015

@author: Alex Lim

Given an input X (which can be modified into a number of different inputs 
[see the file, inputs.py]): X -activates-> Y and X -represses-| Z by the hill
equation. 
"""

import matplotlib.pyplot as plt
import numpy as np
import inputs
from scipy.integrate import odeint
plt.ion()

# Initialize constants, the steady state concentration = B/a
inic = [0,0]    # initial conditions
n = 2.          # hill coefficient
By = 1.         # rate of production of Y
ay = 0.5        # degradation rate of Y
Ky = 1          # activation coefficient of Y
Bz = 1.         # rate of production of Z
az = 0.5        # degradation rate of Z
Kz = 1          # repression coefficient of Z

# Plot dimensions
xMin = 0
xMax = 10
yMin = 0
yMax = 3

itr = 1000      # time iterations
t = np.linspace(xMin, xMax, itr)   # time grid

# Input signal and digital representation
per = 4         # period of input wave
amp = 1         # amplitude of input wave
dis = 1         # displacement of wave (right-shift)
X = inputs.sinInput(t,per,amp,0)
X2 = inputs.squInput(t,per,amp,dis)

# System of DEs and input equation.
def f(ic, t):
        Xi = inputs.sinInput(t,per,amp,0)
        Yi = ic[0]
        Zi = ic[1]    
        # Activator Equation
        f0 = By/(1+np.power(Ky/Xi,n)) - ay*Yi
        # Repressor Equation        
        f1 = Bz/(1+np.power(Xi/Kz,n)) - az*Zi
        return [f0 , f1]

# Solves the system of DEs
soln = odeint(f, inic, t)
Y = soln[:,0]
Z = soln[:,1]

# Plots 
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
# Optional plots (uncomment to have plotted)
# Input Plot of X and digital signal (similar to X when n>>1)
plt.figure()
plt.axis([xMin,xMax,yMin,yMax])
plt.plot(t,X,'b',label = 'input')
plt.plot(t,X2,'b--',label = 'digital')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Input: X concentration')
plt.legend(loc=0)

# Output plot of activator and repressor genes
plt.figure()
plt.axis([xMin,xMax,0,3])
plt.plot(t,Y, 'g', label = 'activator')
plt.plot(t,Z, 'r', label = 'repressor')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)
'''