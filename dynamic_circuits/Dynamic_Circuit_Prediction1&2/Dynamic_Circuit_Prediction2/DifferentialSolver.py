# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 18:19:43 2015

@author: Arinze
"""
import numpy as np

def differentialSolver(func,initc,timeArray,args=()):
    """
    Takes in a function that takes in an array and returns the change 
    in the values as an array, an array of initial conditions, an array of
    time values, and any other arguments necessary for the function.
    
    Returns an array of arrays where the ith index of the returned array 
    corresponds to the values at the ith index of the time array.
    
    This is a back up to odeint that can also be used for comparison. 
    It doesn't use any fancy tricks. It only calculates the linear
    approximation of a value given some set of time values.
    """
    #Convert everything into a numpy array
    results = np.asarray(initc)
    initc = np.asarray(initc)
    theRange = range(1,len(timeArray))
    #Perform the operation once for each time point.
    for i in theRange:
        #get the change in time
        dt = timeArray[i]-timeArray[i-1]
        #use the given function to get the changes in concentrations
        change = np.asarray(func(initc,timeArray[i],*args))
        #update the concentrations for this current time step
        initc = initc + change*dt
        #add the array of current concentrations to the results
        results = np.vstack((results,initc))
    return results