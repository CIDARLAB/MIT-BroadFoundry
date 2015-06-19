# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:36:54 2015

@author: Alex Lim
"""

# cell transcription circuit modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.ion()

#Autoregulators
xMax = 20.
yMax = 12.
B = 1                          # rate of production of X
B0 = 0.05                          # basal production rate
a = 0.1                         # rate of degradation/dilution
K = 5                      # repression threshold
c0 = [K-1,K,K+1,0,K,2*K]  # initial concentration of X
Tr = K/(2*B)                    # response time
n = 50.                         # hill coefficient

# solve the system dy/dt = f(y, t)
def f(ic, t):
        Xi = [0]*3
        Yi = [0]*3
        Xi[0] = ic[0]
        Xi[1] = ic[1]
        Xi[2] = ic[2]
        Yi[0] = ic[3]
        Yi[1] = ic[4]
        Yi[2] = ic[5]
        
        # the model equations
        f0 = B0 + B*(np.power(Xi[0],n)/(np.power(Xi[0],n)+np.power(K,n))) - a*Xi[0]
        f1 = B0 + B*(np.power(Xi[1],n)/(np.power(Xi[1],n)+np.power(K,n))) - a*Xi[1]
        f2 = B0 + B*(np.power(Xi[2],n)/(np.power(Xi[2],n)+np.power(K,n))) - a*Xi[2]
        f3 = B0 + B/(1+(np.power(Yi[0],n)/np.power(K,n))) - a*Yi[0]
        f4 = B0 + B/(1+(np.power(Yi[1],n)/np.power(K,n))) - a*Yi[1]
        f5 = B0 + B/(1+(np.power(Yi[2],n)/np.power(K,n))) - a*Yi[2]
        return [f0,f1,f2,f3,f4,f5]
        
t = np.linspace(0, xMax, 1000) # time grid

# solve the DEs
soln = odeint(f, c0, t)
X1 = soln[:,0]
X2 = soln[:,1]
X3 = soln[:,2]
Y1 = soln[:,3]
Y2 = soln[:,4]
Y3 = soln[:,5]

plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,X1)
plt.plot(t,X2)
plt.plot(t,X3)
plt.axhline(y=K,color='r',ls='--')
plt.title('Autoregulator Activation')
plt.xlabel('Time')
plt.ylabel('Concentration X')
plt.legend(loc=0)

plt.figure()
plt.axis([0, xMax, 0, yMax])
plt.plot(t,Y1)
plt.plot(t,Y2)
plt.plot(t,Y3)
plt.axhline(y=K,color='r',ls='--')
plt.title('Autoregulator Repression')
plt.xlabel('Time')
plt.ylabel('Concentration X')
plt.legend(loc=0)