# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 13:27:25 2015

@author: Alex Lim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

#Input
parts = {'p1':['promoter', 1, 3, 0.1, 0, 0.5, 2, 'p3', 'rep', 'p2', 'rep'], 
         'p2':['promoter', 0, 3, 0.1, 0, 0.5, 2, 'p1', 'rep', 'p3', 'rep'],
         'p3':['promoter', 0, 3, 0.1, 0, 0.5, 2, 'p2', 'rep', 'p1', 'rep']}


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
oTyp = ini[:]
der = ini[:]

i=0
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

# solve the system dx/dt = f(x, t)
def f(init, t):
    Xi = init[0]
    Yi = init[1]
    Zi = init[2]
    der=[Xi,Yi,Zi]
    i=0
    for p in parts:
        if iTyp[i] == 'rep': # Repressor Equation
            X = parts[inp[i]][1]
            if i == 0:
                der[i] = B[i]/(1+np.power(X/K[i],n[i])) - a[i]*der[0]
            elif i == 1:
                der[i] = B[i]/(1+np.power(X/K[i],n[i])) - a[i]*der[1]
            elif i == 2:
                der[i] = B[i]/(1+np.power(X/K[i],n[i])) - a[i]*der[2]
        elif iTyp[i] == 'act': # Activator Equation
            X = parts[inp[i]][1]
            der[i] = B[i]*(np.power(X,n[i])/(np.power(X,n)+np.power(K[i],n[i]))) - a[i]*Y           
        i+=1
    return der

itr = 1000
t = np.linspace(0, 1000, itr)   # time grid

# solve the DEs
soln = odeint(f, ini, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]

# Output: Plot of outputs: activity (activator and repressor)
plt.figure()
plt.plot(t,X, label = 'X')
plt.plot(t,Y, label = 'Y')
plt.plot(t,Z, label = 'Z')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)
