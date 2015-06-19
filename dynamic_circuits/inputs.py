# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:56:25 2015

@author: Alex Lim
Output of not gate IPTG vs activation
Object, object
activating promoter and gene
simulate time scales
production (B) and degradation (a) of generic proteins, rates (# protein/s)
potential transcriptional (assume faster than translation) and translational steps
Promoters models, attributes (inhibition function)
Rbs = factor for the rate, (IPTG vs YFP)-->(t vs YFP @ IPTG = N)

Information
rate Repressor, steady state
value of production rate (Rbs) and degradation rate
Repressor (Pmin, Pmax, K, n)
Gene (alpha(RBS),gamma
design and derivative function for design, generic
defining it
Truth table and inputs

Start with differential equations
concentration (vector)

promoter - not - promoter
then nor
output derivative calculates (only) change of each gene, given current state
treat independently, observe changes

long term goal:
handle circuit with feedback
data structure, dictionary, gene, states and how promoter activates it, where state comes from
external or other gene input, parameter

"""
import numpy as np
from scipy import signal
'''
Sinusoidal Input Function 
per = period of sinusoid 
amp = amplitude of signal
dis = angular displacement to the right
bas = basal production rate
'''
def sinInput(t, per=4, amp=0.8, dis=0, bas=0):
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

'''
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

def squInput(t):
    amp = 3
    per = 2
    return 0.5*amp*(1+signal.square(2*np.pi*t/per))        
    
def sawInput(t):
    amp = 3
    per = 2
    return 0.5*amp*(1+signal.sawtooth(2*np.pi*t/per))
    
    
    
    