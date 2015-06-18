# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:56:25 2015

@author: Alex Lim
"""
import numpy as np

# sinusoidal Input Function (xMax = per*nwaves; yMax = amp+0.1; nwaves = 2)
def sinInput(t):
    per = 4     # period of sinusoid 
    amp = 2     # amplitude of signal
    return 0.5*amp*(1-np.cos(np.pi*t/(0.5*per)))

# Linear Input Function (xMax = 10; yMax = max(xMax*slp,0.5*xMax))  
def linInput(t):
    slp = 1     # slope of line
    intc = 0    # intercept of signal
    return slp * t + intc
    
def stepInput(t, lo, up):
    loB = lo    # lower bound
    upB = up    # upper bound
    step = 2    # magnitude of signal
    if isinstance(t,float):
        if t >= loB and t <= upB:    
            return step
        else:
            return 0
    else:       # For plotting an array x derived from array t
        output = [0]*t.size
        inc=0
        for i in t:
            if i >= loB and i <= upB:
                output[inc] = step
                inc+=1
            else:
                output[inc] = 0
                inc+=1
        return output