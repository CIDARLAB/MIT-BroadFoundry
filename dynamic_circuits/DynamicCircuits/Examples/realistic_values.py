# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:37:25 2015

@author: Alex Lim

mRNA and proteins included in circuit
Values and equations taken and slightly modified from (Elowitz & Leibler, 2000) 
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.ion()

SEC_PER_MIN = 60
MIN_PER_HR = 60
# CONSTANTS
M_HALFLIFE = 2.             # half-life of mRNA transcripts (minutes)
P_HALFLIFE = 10             # half-life of protein (minutes)
Bx = 20.                    # translation efficiency (20 proteins/transcript)
ax = np.log(2)/P_HALFLIFE   # degradation rate of protein (half-life = 10 min)
mBx = 0.4995*SEC_PER_MIN    # induced rate of transcription (30 transcript/min)
mBx0 = 0.0005*SEC_PER_MIN   # repressed rate of transcription (0.03 transcript/min)
amx = np.log(2)/M_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
n = 2                   # hill coefficient
Km = 40                 # repression threshold of promoters


t = np.linspace(xMin, xMax, itr)   # time grid

# solve the system dx/dt = f(x, t)
def f(ic, t):
        Xi = ic[0]
        mX = ic[1]
        
        # Repressor Protein Equations        
        fp0 = Bx*mX - ax*Xi
        # Transcript mRNA Equations
        fm0 = -amx*mX + mBx/(1+np.power(Zi/Km,n)) + mBx0
        fm1 = -amy*mY + mBy/(1+np.power(Xi/Km,n)) + mBy0
        return [fp0, fm0,fm1]
        
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
Xm = soln[:,1]

'''
Reference(s):
Elowitz, M., & Leibler, S. (2000). A synthetic oscillatory network of 
transcriptional regulators. Nature, 403(6767), 335-338. doi:10.1038/35002125
'''