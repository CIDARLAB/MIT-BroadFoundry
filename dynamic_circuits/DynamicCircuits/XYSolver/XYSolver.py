# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
Permutation of interaction graphs between X and Y.
"""

import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.integrate import odeint
plt.ion()


def f(state, t, parameter, vectors):
    '''
    Function passed through odeint() that represents the system of differential 
    equations of the rate of change of each protein.
    '''
    keys = parameter.keys()
    keys.sort()
    initProtein = {}
    
    # Set minimum value to 0.0 if negative and concentration of protein at this state.
    for i in range(len(state)):
        if state[i]<0 :
            state[i] = 0.0
        initProtein[keys[i]] = state[i] 
            
    # Initializes the lists that hold the change in values (deltaValues).
    dinitProtein_dt = []

    # For each protein we want to calculate the change in its levels
    for name in keys:
        # Gets the changes and adds it to the list holding the changes
        dinitProtein_dt.append(getProteinChange(name,initProtein,parameter,vectors))
    return dinitProtein_dt
    
def getFromJsonFile(fileLoc):
    '''
    Reads in a file location and returns the values necessary to generate a graph.
    '''
    myFile = open(fileLoc,'r')
    data = json.load(myFile)
    X_Y_inic = data[0] # a dictionary of initial conditions
    X_Y_params = data[1] # a dictionary containing 2 dictionaries of protein parameters
    X_Y_vectors = data[2:] # a list of lists containing dictionaries of 'to','from', 'effect' vectors
    return X_Y_inic, X_Y_params, X_Y_vectors

# proteinName is a string of the variable that rate of change is being measured
# initProtein is a dictionary of the initial conditions, 
# parameter is a dictionary of 2 dictionaries of protein parameters, 
# vectors is a single list of 4 dictionaries of 'to', 'from' and 'effect'
def getProteinChange(proteinName, initProtein, XYparameter, vectors):
    param = XYparameter[proteinName]

    # Degradation term and basal rate terms
    firstTerm = -param['a']*initProtein[proteinName] + param['B0']
    
    # Regulation term calculation determined by repression or activation hill equations
    name = param['name']
    B = param['B']
    Km = param['Km']
    n = param['n']
    addTerm = 0
    
    if name == 'X':             # Activity of node X
        for vec in vectors:     # Cycles through all vectors finding those that affect X
            if vec['to'] == 'X':
                inputVal = initProtein[vec['from']]
                if vec['effect'] == -1:
                    addTerm += repressor(B,Km,n,inputVal)
                elif vec['effect'] == 1:
                    addTerm += activator(B,Km,n,inputVal)
                else:
                    addTerm += 0
        answer = firstTerm+addTerm
    elif name == 'Y':           #Activity of node Y
        for vec in vectors:
            if vec['to'] == 'Y':
                inputVal = initProtein[vec['from']]
                if vec['effect'] == -1:
                    addTerm += repressor(B,Km,n,inputVal)
                elif vec['effect'] == 1:
                    addTerm += activator(B,Km,n,inputVal)
                else:
                    addTerm += 0
        answer = firstTerm+addTerm
    return answer

#Repressor and Activator Hill Equations
def repressor(B,Km,n,inputVal):
    '''
        Repressor hill equation
    '''
    try:
        return B/(1+(inputVal/Km)**n)
    # If Km is zero
    except ZeroDivisionError:
        try:
            return B*(Km**n)/(inputVal**n+Km**n)
        # If both Km and Input are zero
        except ZeroDivisionError:
            return 0

def activator(B,Km,n,inputVal):
    '''
        Activator hill equation
    '''
    try:
        return B/(1+(Km/inputVal)**n)
    # If Input is zero
    except ZeroDivisionError:
        try:
            return B*(inputVal**n)/(inputVal**n+Km**n)
        # If both Km and Input are zero
        except ZeroDivisionError:
            return 0
            

# Code for plotting configurations of X/Y nodes interactions

# Name of the file
fileLoc = 'masterXYInputFile4.json'
xMax = 200
itr = 1000
t = np.linspace(0, xMax, itr)
X_Y_initCond, X_Y_params, X_Y_vectors = getFromJsonFile(fileLoc)
initc = []
keys = X_Y_initCond.keys()
keys.sort()
for key in keys:
    initc.append(X_Y_initCond[key])
    
for i in range(len(X_Y_vectors)):
    soln_f = odeint(f, initc, t, args = (X_Y_params, X_Y_vectors[i]))
    X = soln_f[:,0]
    Y = soln_f[:,1]
    print X_Y_vectors[i]
    plt.figure()
    #plt.axis([xMin, xMax, yMin, yMax])
    plt.plot(t,X, 'r--', label = 'X')
    plt.plot(t,Y, 'b--', label = 'Y')
    plt.xlabel('Time (min)')
    plt.ylabel('Proteins per cell')
    plt.title('Example protein data')
    plt.legend(loc=0)
