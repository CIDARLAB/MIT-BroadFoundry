# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:37:25 2015

@author: Alex Lim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

xMin = 0                # plot dimensions
xMax = 1000
yMin = 0
yMax = 5000
itr = 1000              # number of time iterations
inic = [20,0,0,0,0,0]   # initial system conditions
Bx = 20.                # translation efficiency (protein/transcript)
By = 20.
Bz = 20.
ax = np.log(2)/10       # degradation rate of protein (half-life = 10 min)
ay = np.log(2)/10
az = np.log(2)/10
mBx = 0.4995*60         # induced rate of transcription (transcript/min)
mBy = 0.4995*60
mBz = 0.4995*60
mBx0 = 0.0005*60        # repressed rate of transcription (transcript/min)
mBy0 = 0.0005*60
mBz0 = 0.0005*60
amx = np.log(2)/2.       # degradation rate of mRNA (half-life = 2 min)
amy = np.log(2)/2.
amz = np.log(2)/2.


n = 2                   # hill coefficient
Km = 40                 # repression threshold of X

t = np.linspace(xMin, xMax, itr)   # time grid

# solve the system dx/dt = f(x, t)
def f(ic, t):
        Xi = ic[0]        
        Yi = ic[1]        
        Zi = ic[2]
        mX = ic[3]
        mY = ic[4]            
        mZ = ic[5]        
        # Repressor Protein Equations        
        fp0 = Bx*mX - ax*Xi
        fp1 = By*mY - ay*Yi
        fp2 = Bz*mZ - az*Zi
        # Transcript mRNA Equations
        fm0 = -amx*mX + mBx/(1+np.power(Zi/Km,n)) + mBx0
        fm1 = -amy*mY + mBy/(1+np.power(Xi/Km,n)) + mBy0
        fm2 = -amz*mZ + mBz/(1+np.power(Yi/Km,n)) + mBz0
        return [fp0, fp1, fp2, fm0, fm1, fm2]
                
# solve the DEs
soln = odeint(f, inic, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]

plt.figure()
plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,X, label = 'X')
plt.plot(t,Y, label = 'Y')
plt.plot(t,Z, label = 'Z')
plt.xlabel('Time (min)')
plt.ylabel('Protein per cell')
plt.title('Example Data')
plt.legend(loc=0)