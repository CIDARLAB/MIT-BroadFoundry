# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:33:26 2015

@author: Arinze
"""

import matplotlib.pyplot as plt
import inputs
from scipy.integrate import odeint
import numpy as np

xMin = 0
#yMin = 0
xMax = 2000
#yMax = ?
itr = 500      # time iterations
t = np.linspace(xMin, xMax, itr)   # time grid

mXi = 10
mYi = 10

Xi = 0
Yi = 0

initc = [mXi, Xi, mYi, Yi]

Mx_HALFLIFE = 2.             # half-life of mRNA transcripts
Px_HALFLIFE = 10             # half-life of protein
My_HALFLIFE = 2.             # half-life of mRNA transcripts
Py_HALFLIFE = 10             # half-life of protein

amx = np.log(2)/Mx_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
mBx1 = 30    # induced rate of transcription (30 transcript/min) by IPTG
mBx2 = 30    # induced rate of transcription (30 transcript/min) by IPTG
amy = np.log(2)/My_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
mBy = 30    # induced rate of transcription (30 transcript/min)

nx1 = 2                   # hill coefficient for protein IPTG
Kmx1 = 40                 # repression threshold of promoter for IPTG, monomers per cell
nx2 = 2                   # hill coefficient for protein IPTG2
Kmx2 = 40                 # repression threshold of promoter for IPTG2, monomers per cell
ny = 2                   # hill coefficient
Kmy = 30000                 # repression threshold of promoters, monomers per cell

#IPTG_low = inputs.linInput(t, 0, 0.1)

Bx = 20.                    # translation efficiency (20 proteins/transcript)
ax = np.log(2)/Px_HALFLIFE   # degradation rate of protein (half-life = 10 min)
By = 20.                    # translation efficiency (20 proteins/transcript)
ay = np.log(2)/Px_HALFLIFE   # degradation rate of protein (half-life = 10 min)

var1 = [inputs.squInput,1000,100,0,0]
var2 = [inputs.squInput,500,100,0,0]
#f1 = inputs.squInput
#f2 = inputs.squInput

# System of DEs and input equation.
def f(state, t):
        #concentration of mRNA at this state
        mX = state[0]
        Xi = state[1]        
        mY = state[2]
        Yi = state[3]
        
        #IPTG input
        IPTGi = var1[0](t,var1[1],var1[2],var1[3],var1[4])
        
        #IPTG2 input
        IPTGi2 = var2[0](t,var2[1],var2[2],var2[3],var2[4])
        
        #mRNA production
        try:
            dmX_dt = -amx*mX + mBx1*(Kmx1**nx1)/(IPTGi**nx1+Kmx1**nx1) + mBx2*(Kmx2**nx2)/(IPTGi2**nx2+Kmx2**nx2)
            dmY_dt = -amy*mY + mBy*(Kmy**ny)/(Xi**ny+Kmy**ny)
        except OverflowError:
            dmX_dt = -amx*mX + mBx1/(1+(IPTGi/Kmx1)**nx1) + mBx2/(1+(IPTGi2/Kmx2)**nx2)
            dmY_dt = -amy*mY + mBy/(1+(Xi/Kmy)**ny)
        #protein production
        dX_dt = Bx*mX - ax*Xi
        dY_dt = By*mY - ay*Yi
        
        return [dmX_dt , dX_dt , dmY_dt , dY_dt]
        
        

soln_f = odeint(f, initc, t)
#soln_g = odeint(g, initc, t)

Xm_f = soln_f[:,0]
Xp_f = soln_f[:,1]
Ym_f = soln_f[:,2]
Yp_f = soln_f[:,3]

#Xm_g = soln_g[:,0]
#Xp_g = soln_g[:,1]

#Visualizing IPTG
IPTG = []
IPTG2 = []
for timepos in t:
    IPTGi = var1[0](timepos,var1[1],var1[2],var1[3],var1[4])
    IPTGi2 = var2[0](timepos,var2[1],var2[2],var2[3],var2[4])
    IPTG.append(IPTGi)
    IPTG2.append(IPTGi2)
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.axis([0, 2000, 0, 120])
plt.plot(t,IPTG, 'r--', label = 'IPTG')
plt.plot(t,IPTG2, 'r--', label = 'IPTG2')
plt.xlabel('Time (min)')
plt.ylabel('IPTG')
plt.title('IPTG present')
plt.legend(loc=0)

#Visualizing mRNA
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,Xm_f, 'g--', label = 'Xm_f')
plt.plot(t,Ym_f, 'm--', label = 'Ym_f')
plt.xlabel('Time (min)')
plt.ylabel('Transcripts per cell')
plt.title('Example mRNA Data')
plt.legend(loc=0)

#Visualizing protein
plt.figure()
#plt.axis([xMin, xMax, 0, 500])
plt.plot(t,Xp_f, 'g--', label = 'X_f')
plt.plot(t,Yp_f, 'm--', label = 'Y_f')
plt.xlabel('Time (min)')
plt.ylabel('Proteins per cell')
plt.title('Example Protein Data')
plt.legend(loc=0)


