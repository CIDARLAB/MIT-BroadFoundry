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
xMax = 2000


mXi = 10

Xi = 0

initc = [mXi, Xi]

Mx_HALFLIFE = 2.             # half-life of mRNA transcripts
Px_HALFLIFE = 10             # half-life of protein

amx = np.log(2)/Mx_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
mBx = 30    # induced rate of transcription (30 transcript/min)


n = 2                   # hill coefficient
Km = 40                 # repression threshold of promoters, monomers per cell

itr = 500      # time iterations
t = np.linspace(xMin, xMax, itr)   # time grid

#IPTG_low = inputs.linInput(t, 0, 0.1)

Bx = 20.                    # translation efficiency (20 proteins/transcript)
ax = np.log(2)/Px_HALFLIFE   # degradation rate of protein (half-life = 10 min)


# System of DEs and input equation.
def f(state, t):
        #concentration of mRNA at this state
        mX = state[0]
        Xi = state[1]
        
        #IPTG input
        #Calculates the value of IPTG at that time.
        #IPTGi = inputs.linInput(t, 0, 0.5)
        #IPTGi = inputs.stepFunction(t, 500, 100, 0.01)
        IPTGi = inputs.sinInput(t, 500, 100, 0.01)
        #IPTGi = inputs.squInput(t, per=500, amp=100, x0=0, bas=0)
        
        #mRNA production
        dmX_dt = -amx*mX + mBx*(Km**n)/(IPTGi**n+Km**n)

        
        #protein production
        dX_dt = Bx*mX - ax*Xi    
        
        return [dmX_dt , dX_dt]
        
        

soln_f = odeint(f, initc, t)
#soln_g = odeint(g, initc, t)

Xm_f = soln_f[:,0]
Xp_f = soln_f[:,1]

#Xm_g = soln_g[:,0]
#Xp_g = soln_g[:,1]



#Visualizing mRNA
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,Xm_f, 'r--', label = 'Xm_f')
plt.plot(t,Xp_f/250, 'r-', label = 'X_f')
plt.xlabel('Time (min)')
plt.ylabel('Transcripts per cell')
plt.title('Example mRNA Data')
plt.legend(loc=0)

#Visualizing protein
plt.figure()
#plt.axis([xMin, xMax, 0, 500])
plt.plot(t,Xp_f, 'r-', label = 'X_f')
plt.xlabel('Time (min)')
plt.ylabel('Proteins per cell')
plt.title('Example Protein Data')
plt.legend(loc=0)

#Visualizing IPTG
IPTG = []
for timepos in t:
    #IPTGi = inputs.linInput(timepos, 0, 0.5)
    #IPTGi = inputs.stepFunction(timepos, 500, 100, 0.01)
    IPTGi = inputs.sinInput(timepos, 500, 100, 0.01)
    #IPTGi = inputs.squInput(timepos, per=500, amp=100, x0=0, bas=0)
    IPTG.append(IPTGi)    
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,IPTG, 'b--', label = 'IPTG')
plt.xlabel('Time (min)')
plt.ylabel('IPTG')
plt.title('IPTG present')
plt.legend(loc=0)

