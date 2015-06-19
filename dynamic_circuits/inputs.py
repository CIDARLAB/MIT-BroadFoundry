# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:56:25 2015

@author: Alex Lim
"""
import numpy as np
from scipy import signal

# sinusoidal Input Function (xMax = per*nwaves; yMax = amp+0.1; nwaves = 2)
def sinInput(t, per=4, amp=2, dis=0, bas=0):
    per = 4     # period of sinusoid 
    amp = 2     # amplitude of signal
    dis = 0
    return 0.5*amp*(1-np.cos(np.pi*(t-dis)/(0.5*per)))

# Linear Input Function (xMax = 10; yMax = max(xMax*slp,0.5*xMax))  
def linInput(t):
    slp = 1     # slope of line
    intc = 0    # intercept of signal
    return slp * t + intc
    
def stepInput(t, lo, up): # variation 1
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

def squInput(t):
    amp = 3
    per = 2
    return 0.5*amp*(1+signal.square(2*np.pi*t/per))        
    
def sawInput(t):
    amp = 3
    per = 2
    return 0.5*amp*(1+signal.sawtooth(2*np.pi*t/per))
    
    
    
    