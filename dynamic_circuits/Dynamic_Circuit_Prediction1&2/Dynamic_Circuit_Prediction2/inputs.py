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
    if type(t)==list:
        t = np.asarray(t)
    return bas + 0.5*amp*(1-np.cos(np.pi*(t-dis)/(0.5*per)))
    
'''
Linear Input Function 
slp = slope of line
intc = intercept of signal
'''
def linInput(t, slp=0, intc=1):
    if type(t)==list:
        t = np.asarray(t)
    return slp * t + intc

'''
Single Step Input Function
x0 = starting x-position of step signal
step = magnitude of step signal
bas = base of step signal
'''
def stepFunction(t,x0=3,step=2, bas=0):
    if isinstance(t,float):
        if t >= x0:    
            return bas+step
        else:
            return bas
    else:       # For plotting an array x derived from array t
        output = [0]*t.size
        inc=0
        for i in t:
            if i >= x0:
                output[inc] = bas+step
                inc+=1
            else:
                output[inc] = bas
                inc+=1
        return output
    
'''
Single Step Up and Down Input Function
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
def squInput(t, per=2, amp=1, x0=0, bas=0):
    if type(t)==list:
        t = np.asarray(t)
    return bas + 0.5*amp*(1+signal.square(2.0*np.pi*(t-x0)/per))

'''
Sawtooth Wave Input Function
amp = Amplitude of wave
per = Period of wave
bas = base of wave
'''    
def sawInput(t, amp=3, per=2, bas=0):
    if type(t)==list:
        t = np.asarray(t)
    return bas + 0.5*amp*(1+signal.sawtooth(2*np.pi*t/per))

'''
Square wave where you can specify when it is high or low within one period
maxVal = value at a high
minVal = value at a min
per = period of wave
tv = specification of when it should be high. String of 1s and 0s
'''
def specInput(t, maxVal=1, minVal=0, per=1.0, tv="01"):
    if type(t)==list:
        t = np.asarray(t)
    numSegments = len(tv)
    increment = float(per)/numSegments
    t = t%per
    
    if type(t)!=np.ndarray:
        for segment in range(1,numSegments+1):
            if t<increment*segment:
                segment-=1
                break
        
        value = tv[segment] 
        if value=="1":
            return maxVal
        elif value=="0":
            return minVal
        else:
            print "Error truth value not properly formatted"
            return

    output = [0]*t.size
    for i in range(t.size):
        tval = t[i]
    
        for segment in range(1,numSegments+1):
            if tval<increment*segment:
                segment-=1
                break
        value = tv[segment] 
        if value=="1":
            output[i] = maxVal
        elif value=="0":
            output[i] = minVal
        else:
            print "Error truth value not properly formatted"
            return
    return output
    