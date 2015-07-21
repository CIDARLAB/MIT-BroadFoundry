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

mZi = 10

Zi = 0

initc = [mZi, Zi]

Mz_HALFLIFE = 2.             # half-life of mRNA transcripts
Pz_HALFLIFE = 10             # half-life of protein

amz = np.log(2)/Mz_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
mBz1 = 30    # induced rate of transcription (30 transcript/min) by IPTG
mBz2 = 30    # induced rate of transcription (30 transcript/min) by IPTG2

nz1 = 2                   # hill coefficient for protein IPTG
Kmz1 = 40                 # repression threshold of promoter for IPTG, monomers per cell
nz2 = 2                   # hill coefficient for protein IPTG2
Kmz2 = 40                 # repression threshold of promoter for IPTG2, monomers per cell

itr = 500      # time iterations
t = np.linspace(xMin, xMax, itr)   # time grid

#IPTG_low = inputs.linInput(t, 0, 0.1)

Bz = 20.                    # translation efficiency (20 proteins/transcript)
az = np.log(2)/Pz_HALFLIFE   # degradation rate of protein (half-life = 10 min)

var1 = [inputs.squInput,1000,100,0,0]
var2 = [inputs.squInput,500,100,0,0]
#f1 = inputs.squInput
#f2 = inputs.squInput

# System of DEs and input equation.
def f(state, t):
        #concentration of mRNA at this state
        mZ = state[0]
        Zi = state[1]        
        
        #IPTG input
        IPTGi = var1[0](t,var1[1],var1[2],var1[3],var1[4])
        
        #IPTG2 input
        IPTGi2 = var2[0](t,var2[1],var2[2],var2[3],var2[4])
        
        #mRNA production
        try:
            dmZ_dt = -amz*mZ + mBz1*(Kmz1**nz1)/(IPTGi**nz1+Kmz1**nz1) + mBz2*(Kmz2**nz2)/(IPTGi2**nz2+Kmz2**nz2)
        except OverflowError:
            dmZ_dt = -amz*mZ + mBz1/(1+(IPTGi/Kmz1)**nz1) + mBz2/(1+(IPTGi2/Kmz2)**nz2)
        #protein production
        dZ_dt = Bz*mZ - az*Zi
        
        return [dmZ_dt , dZ_dt]
        
        

soln_f = odeint(f, initc, t)
#soln_g = odeint(g, initc, t)

Zm_f = soln_f[:,0]
Zp_f = soln_f[:,1]

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
plt.plot(t,IPTG2, 'b--', label = 'IPTG2')
plt.xlabel('Time (min)')
plt.ylabel('IPTG')
plt.title('IPTG present')
plt.legend(loc=0)

#Visualizing mRNA
plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,Zm_f, 'g--', label = 'Zm_f')
plt.xlabel('Time (min)')
plt.ylabel('Transcripts per cell')
plt.title('Example mRNA Data')
plt.legend(loc=0)

#Visualizing protein
plt.figure()
#plt.axis([xMin, xMax, 0, 500])
plt.plot(t,Zp_f, 'g--', label = 'Z_f')
plt.xlabel('Time (min)')
plt.ylabel('Proteins per cell')
plt.title('Example Protein Data')
plt.legend(loc=0)


