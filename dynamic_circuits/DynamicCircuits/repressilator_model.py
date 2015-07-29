# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:48:34 2015

@author: Alex Lim
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

xMin = 0
xMax = 300
yMin = 0
yMax = 3
itr = 1000          # time iterations
inic = [2,0,0]      # initial conditions
Bx = 10              # rate of production of X
ax = 0.1            # degradation rate of X
Bx0 = 0             # basal production rate of X
By = 10              # rate of production of Y
ay = 0.1            # degradation rate of Y
By0 = 0             # basal production rate of Y
Bz = 10              # rate of production of Z
az = 0.1            # degradation rate of Z
Bz0 = 0             # basal production rate of Z
n = 2               # hill coefficient
Kx = 1            # repression coefficient of X
Ky = 1            # repression coefficient of Y
Kz = 1            # repression coefficient of Z

t = np.linspace(xMin, xMax, itr)   # time grid

# solve the system dx/dt = f(x, t)
def f(ic, t):
        Xi = ic[0]
        Yi = ic[1]    
        Zi = ic[2]        
        # Repressor Equations
        f0 = Bx/(1+np.power(Zi/Kx,n)) - ax*Xi + Bx0
        f1 = By/(1+np.power(Xi/Ky,n)) - ay*Yi + By0       
        f2 = Bz/(1+np.power(Yi/Kz,n)) - az*Zi + Bz0
        return [f0, f1, f2]
                
# solve the DEs
soln = odeint(f, inic, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]

plt.figure()
plt.plot(t,X, label = 'X')
plt.plot(t,Y, label = 'Y')
plt.plot(t,Z, label = 'Z')
plt.xlabel('time')
plt.ylabel('concentration')
plt.title('Output: Activity')
plt.legend(loc=0)





