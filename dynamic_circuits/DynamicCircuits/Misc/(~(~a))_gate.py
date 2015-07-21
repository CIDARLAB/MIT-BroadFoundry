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
mBx = 30    # induced rate of transcription (30 transcript/min)
amy = np.log(2)/My_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
mBy = 30    # induced rate of transcription (30 transcript/min)

nx = 999                   # hill coefficient
Kmx = 40                 # repression threshold of promoters, monomers per cell
ny = 999                  # hill coefficient
Kmy = 40                 # repression threshold of promoters, monomers per cell

itr = 500      # time iterations
t = np.linspace(xMin, xMax, itr)   # time grid

#IPTG_low = inputs.linInput(t, 0, 0.1)

Bx = 20.                    # translation efficiency (20 proteins/transcript)
ax = np.log(2)/Px_HALFLIFE   # degradation rate of protein (half-life = 10 min)
By = 20.                    # translation efficiency (20 proteins/transcript)
ay = np.log(2)/Py_HALFLIFE   # degradation rate of protein (half-life = 10 min)

# System of DEs and input equation.
def f(state, t):
        #concentration of mRNA at this state
        mX = state[0]
        Xi = state[1]
        mY = state[2]
        Yi = state[3]
        
        #IPTG input
        #Calculates the value of IPTG at that time.
        #IPTGi = inputs.linInput(t, 0, 0.5)
        #IPTGi = inputs.stepFunction(t, 500, 1000, 0.01)
        IPTGi = inputs.sinInput(t, 500, 1000, 0.01)
        #IPTGi = inputs.squInput(t, per=500, amp=100, x0=0, bas=0)
        
        #mRNA production
        try:
            dmX_dt = -amx*mX + mBx*(Kmx**nx)/(IPTGi**nx+Kmx**nx)
            dmY_dt = -amy*mY + mBy*(Kmy**ny)/(Xi**ny+Kmy**ny)
        except OverflowError:
            dmX_dt = -amx*mX + mBx/(1+(IPTGi/Kmx)**nx)
            dmY_dt = -amy*mY + mBy/(1+(Xi/Kmy)**ny)
        #protein production
        dX_dt = Bx*mX - ax*Xi
        dY_dt = By*mY - ay*Yi
        
        return [dmX_dt , dX_dt, dmY_dt , dY_dt]
        
        

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
for timepos in t:
    #IPTGi = inputs.linInput(timepos, 0, 0.5)
    #IPTGi = inputs.stepFunction(timepos, 500, 1000, 0.01)
    IPTGi = inputs.sinInput(timepos, 500, 1000, 0.01)
    #IPTGi = inputs.squInput(timepos, per=500, amp=100, x0=0, bas=0)
    IPTG.append(IPTGi)    
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,IPTG, 'g--', label = 'IPTG')
plt.xlabel('Time (min)')
plt.ylabel('IPTG')
plt.title('IPTG present')
plt.legend(loc=0)

#Visualizing mRNA
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,Xm_f, 'r--', label = 'Xm_f')
plt.plot(t,Ym_f, 'b--', label = 'Ym_f')
plt.xlabel('Time (min)')
plt.ylabel('Transcripts per cell')
plt.title('Example mRNA Data')
plt.legend(loc=0)

#Visualizing protein
plt.figure()
#plt.axis([xMin, xMax, 0, 500])
plt.plot(t,Xp_f, 'r-', label = 'X_f')
plt.plot(t,Yp_f, 'b--', label = 'Y_f')
plt.xlabel('Time (min)')
plt.ylabel('Proteins per cell')
plt.title('Example Protein Data')
plt.legend(loc=0)


