# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:16:34 2015

@author: Alex Lim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

#Positive Feed-Forward Loop with And Gate
xMax = 20.
yMax = 1
# B0 = 0.                     # basal rate of production 
Bxi = 2.                      # rate of production of X
Byi = 2.                      # rate of production of Y
Bzi = 2.                      # rate of production of Z
a = 0.1                       # rate of degradation/dilution (same for X,Y & Z)
c0 = [0,0,0]                  # initial concentrations [X,Y,Z]
Xst=Bxi/a                     # steady state of X
Yst=Byi/a                     # steady state of Y
Zst=Bzi/a                     # steady state of Z
Kxy = 5.                      # activation threshold of Y
Kyz = 5.                      # activation threshold of Z
xTr = np.log(2)/a             # response time
# c0[0] < Kxy
yTres = max(np.log((Xst-c0[0])/(Xst-Kxy))/a,0)     # Y activation time
zTres = max(np.log((Yst-c0[1])/(Yst-Kyz))/a,0)     # Z activation time
Bx = Bxi
itr = 1000

# solve the system dx/dt = f(x, t)
def f(x, t):
        Xi = x[0]
        Yi = x[1]
        Zi = x[2]
        # the model equations
        f0 = Bx - a*Xi
        f1 = By - a*Yi
        f2 = Bz - a*Zi
        return [f0,f1,f2]
        
tx = np.linspace(0, xMax, itr)   # time grid

if c0[1] != 0 and yTres != 0:
    ty = np.linspace(0, yTres, itr)
    By = 0
    sol = odeint(f, c0,ty)
    Y = sol[:,1]
    c0[1]=Y[itr-1]
By = Byi
ty2 = np.linspace(yTres, xMax, itr)

if c0[2] !=0 and zTres != 0:
    if yTres != 0:    # initial decay
        tz = np.linspace(0, yTres, itr)
        Bz = 0
        sol = odeint(f,c0,tz)
        Z = sol[:,2]
        c0[2]=Z[itr-1]
    tz2 = np.linspace(yTres, yTres+zTres, itr)
    Bz = 0
    sol = odeint(f,c0,tz2)
    Z2 = sol[:,2]        
    c0[2]=Z2[itr-1]
Bz = Bzi
  
tz3 = np.linspace(yTres+zTres, xMax, itr)

# solve the DEs
xSoln = odeint(f, c0, tx)
ySoln = odeint(f, c0, ty2)
zSoln = odeint(f, c0, tz3)
X = xSoln[:,0]
Y2 = ySoln[:,1]
Z3 = zSoln[:,2]

heightX = (c0[0]+(Xst-c0[0])*(1-np.power((np.e),-1*a*yTres)))/(Xst*yMax)
heightY = (c0[1]+(Yst-c0[1])*(1-np.power((np.e),-1*a*zTres)))/(Yst*yMax)

# Three subplots sharing both x/y axes
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

# X
ax1.plot(tx, X/Xst)
ax1.set_xlim([0, xMax])
ax1.set_ylim([0, yMax])
ax1.axvline(x=yTres, ymax=heightX, color='r', ls='--')
ax1.set_ylabel('[X]/[X]st')
ax1.set_title('Positive Feed-Forward Loop (X -> Y, X & Y -> Z)')

# Y
if c0[1] != 0 and yTres != 0:
    ax2.plot(ty,Y/Yst)
ax2.plot(ty2, Y2/Yst)
ax2.axvline(x=yTres, color='r', ls='--')
ax2.axvline(x=yTres+zTres, ymax=heightY, color='r', ls='--')
ax2.set_ylabel('[Y]/[Y]st')

# Z
if c0[2] != 0 and zTres != 0:
    if yTres != 0:
        ax3.plot(tz, Z/Zst)
    ax3.plot(tz2, Z2/Zst)

ax3.plot(tz3, Z3/Zst)
ax3.axvline(x=yTres, color='r', ls='--')
ax3.axvline(x=yTres+zTres, color='r', ls='--')
ax3.set_ylabel('[Z]/[Z]st')

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0)
plt.setp([i.get_xticklabels() for i in f.axes[:-1]], visible=False)
plt.xlabel('time')

