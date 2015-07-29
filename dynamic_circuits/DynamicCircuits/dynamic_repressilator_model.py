# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:37:25 2015

@author: Alex Lim

Numbers and calculations taken and slightly modified from 
(Elowitz & Leibler, 2000) 
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

# CONSTANTS
M_HALFLIFE = 2.             # half-life of mRNA transcripts
P_HALFLIFE = 10             # half-life of protein
SEC_PER_MIN = 60
MIN_PER_HR = 60
 
xMin = 0                    # plot dimensions
xMax = 500
yMin = 0
itr = 1000                  # number of time iterations
inic = [0,30,0,0,0,0]       # initial system conditions
Bx = 20.                    # translation efficiency (20 proteins/transcript)
By = 20.
Bz = 20.
ax = np.log(2)/P_HALFLIFE   # degradation rate of protein (half-life = 10 min)
ay = np.log(2)/P_HALFLIFE
az = np.log(2)/P_HALFLIFE
mBx = 0.4995*SEC_PER_MIN    # induced rate of transcription (30 transcript/min)
mBy = 0.4995*SEC_PER_MIN
mBz = 0.4995*SEC_PER_MIN
mBx0 = 0.0005*SEC_PER_MIN   # repressed rate of transcription (0.03 transcript/min)
mBy0 = 0.0005*SEC_PER_MIN
mBz0 = 0.0005*SEC_PER_MIN
amx = np.log(2)/M_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
amy = np.log(2)/M_HALFLIFE
amz = np.log(2)/M_HALFLIFE

n = 2                   # hill coefficient
Km = 40                 # repression threshold of promoters

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
        
'''
The concentration of protein depends on the concentration of mRNA transcripts 
and degradation of the protein itself.
Concentration of mRNA (X,Y,Z) is repressed by a specific protein (Z,X,Y) 
according to the Hill equation  
Assumption of continuous values rather than discrete values.
'''     
                
# solve the DEs
soln = odeint(f, inic, t)
X = soln[:,0]
Y = soln[:,1]
Z = soln[:,2]
Xm = soln[:,3]
Ym = soln[:,4]
Zm = soln[:,5]

yMax = max([1.1*max(X),1.1*max(Y),1.1*max(Z)]) # plot max Y height for protein

plt.figure()
plt.plot(t,Xm, 'c--', label = 'Xm')
plt.plot(t,Ym, 'y--', label = 'Ym')
plt.plot(t,Zm, 'm--', label = 'Zm')
plt.xlabel('Time (min)')
plt.ylabel('Transcripts per cell')
plt.title('Example mRNA Data')
plt.legend(loc=0)

plt.figure()
plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,X, 'b-', label = 'X')
plt.plot(t,Y, 'g-', label = 'Y')
plt.plot(t,Z, 'r-', label = 'Z')
plt.xlabel('Time (min)')
plt.ylabel('Proteins per cell')
plt.title('Example Protein Data')
plt.legend(loc=0)

# Plot of both
fig, ax1 = plt.subplots()
ax1.plot(t,X, 'b-', label = 'X')
ax1.plot(t,Y, 'g-', label = 'Y')
ax1.plot(t,Z, 'r-', label = 'Z')
ax1.set_xlabel('Time (min)')
ax1.set_ylabel('Proteins per cell')

ax2 = ax1.twinx()
plt.plot(t,Xm, 'c--', label = 'Xm')
plt.plot(t,Ym, 'y--', label = 'Ym')
plt.plot(t,Zm, 'm--', label = 'Zm')
plt.ylabel('Transcripts per cell')


'''
Reference(s):
Elowitz, M., & Leibler, S. (2000). A synthetic oscillatory network of 
transcriptional regulators. Nature, 403(6767), 335-338. doi:10.1038/35002125
'''