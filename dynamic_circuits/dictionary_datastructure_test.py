# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:37:25 2015

@author: Alex Lim
"""

'''
STILL IN THE WORKINGS...
reading paper
code copied from repressilator_model.py
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

'''
Dictionary of parts represented by lists containing characteristics
{k:[name, ic, B, a, B0, inputs, iType, outputs, oType]}  
Use [if s[4] is None:] for boolean
'''
parts = {'p1':['promoter', 0, 5, 0.1, 0, 0.5, 2, 'p3', 'rep', 'p2', 'rep'], 
         'p2':['promoter', 0, 4, 0.1, 0, 0.5, 2, 'p1', 'rep', 'p3', 'rep'],
         'p3':['promoter', 0, 3, 0.1, 0, 0.5, 2, 'p2', 'rep', 'p1', 'rep']}



# solve the system dx/dt = f(x, t)
def f(init, t):
        size = len(parts)
        ini = [0]*size
        B = ini[:]
        a = ini[:]
        B0 = ini[:]
        K = ini[:]
        n = ini[:]
        inp = ini[:]
        iTyp = ini[:]
        outp = ini[:]
        oTyp = inic[:]
        der = inic[:]
        
        i = 0
        # Creates parameters of DEs        
        for p in parts:
            if(parts[p][0] == 'promoter'):
                ini[i]= parts[p][1]
                B[i] = parts[p][2]
                a[i] = parts[p][3]
                B0[i] = parts[p][4]
                K[i] = parts[p][5]
                n[i] = parts[p][6]
                inp[i] = parts[p][7]
                iTyp[i] = parts[p][8]
                outp[i] = parts[p][9]
                oTyp[i] = parts[p][-1]
            i += 1
        
        i = 0
        # Differential Equations for Repressor and Activator
        for p in parts:
            if iTyp[i] == 'rep':
                der[i] = B[i]/(1+np.power(parts[inp[i]][1]/K[i],n[i])) - a[i]*ini[i] + B0[i]
            elif iTyp[i] == 'act':
                der[i] = B[i]*np.power(parts[inp[i]][1], n[i])/(np.power(parts[inp[i]][1],n[i])+np.power(K[i],n[i])) - a[i]*ini[i] + B0[i]
            i += 1
        print der
        return der

xMin = 0
xMax = 300
yMin = 0
yMax = 3
itr = 1000          # time iterations
inic = [0,0,0]      # initial conditions
n = 2               # hill coefficient
Kx = 0.5            # repression coefficient of X
Ky = 0.5            # repression coefficient of Y
Kz = 0.5            # repression coefficient of Z

t = np.linspace(xMin, xMax, itr)   # time grid

                
# solve the DEs
soln = odeint(f, inic, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]

plt.figure()
plt.plot(t,X, label = 'p1')
plt.plot(t,Y, label = 'p2')
plt.plot(t,Z, label = 'p3')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)