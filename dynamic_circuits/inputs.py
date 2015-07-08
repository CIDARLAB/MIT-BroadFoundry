# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:56:25 2015

@author: Alex Lim
"""

import numpy as np
from scipy import signal
'''
Sinusoidal Input Function 
per = period of sinusoid 
amp = amplitude of signal
dis = angular displacement to the right
bas = base of sinusoid wave
'''
def sinInput(t, per=4, amp=1, dis=0, bas=0):
    return bas + 0.5*amp*(1-np.cos(np.pi*(t-dis)/(0.5*per)))
    
'''
Linear Input Function 
slp = slope of line
intc = intercept of signal
'''
def linInput(t, slp=0, intc=1):
    return slp * t + intc
    
'''
Single Step Input Function
x0 = starting x-position of step signal
dis = length of step signal
step = magnitude of step signal
bas = base of step signal
'''
def stepInput(t, x0=3, dis=2, step=2, bas=0):
    if isinstance(t,float):
        if t >= x0 and t <= x0+dis:    
            return bas+step
        else:
            return bas
    else:       # For plotting an array x derived from array t
        output = [0]*t.size
        inc=0
        for i in t:
            if i >= x0 and i <= x0+dis:
                output[inc] = bas+step
                inc+=1
            else:
                output[inc] = bas
                inc+=1
        return output

'''
Square Wave Input Function
amp = Amplitude of square wave
per = Period of square wave
bas = base of wave
'''
def squInput(t, amp=3, per=2, bas=0):
    return bas + 0.5*amp*(1+signal.square(2*np.pi*t/per))        

'''
Sawtooth Wave Input Function
amp = Amplitude of wave
per = Period of wave
bas = base of wave
'''    
def sawInput(t, amp=3, per=2, bas=0):
    return bas + 0.5*amp*(1+signal.sawtooth(2*np.pi*t/per))
    
    
    
    