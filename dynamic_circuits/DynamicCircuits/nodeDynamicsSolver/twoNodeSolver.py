# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
Permutation of interaction graphs between X and Y nodes.

NOTE: 
    json array == [python list] ; 
    json objects == {python dictionary}
"""

import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.integrate import odeint
plt.ion()

# fileLoc: directory of the json file.
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

# state: the initial condition at which the differential would be determined
# t: the list of time values
# parameter: a dictionary of gene parameters where the keys are the gene names
#   -For example, 'X' would be a key to a dictionary of parameters
# vectors: an array of objects with each specifying a regulation motif 
def f(state, t, parameter, vectors):
    '''
    Function passed through odeint() that represents the system of differential 
    equations of the rate of change of each protein.
    '''
    keys = parameter.keys()
    keys.sort()
    initProtein = {}
    
    # By convention, the first index of state corresponds to the first key in the initProtein
    # Set minimum value to 0.0 if the concentration of protein at this state is negative.
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

# proteinName: a string of the variable that rate of change is being measured
# initProtein: a dictionary of the initial conditions, 
# parameter: a dictionary of 2 dictionaries of protein parameters, 
# vectors: a single array of objects representing regulation motifs.
def getProteinChange(proteinName, initProtein, XYparameter, vectors):
    param = XYparameter[proteinName]

    a = param['a']
    B0 = param['B0']
    # Degradation term and basal rate terms
    firstTerm = -a*initProtein[proteinName] + B0
    # Regulation term calculation determined by repression or activation hill equations
    addTerm = 0   
    
    for vec in vectors:     
        # Cycles through all vectors finding those that affect the node
        if vec['to'] == proteinName:
            # Sets up terms for the calculation of the hill equation                                
            inputVal = initProtein[vec['from']]
            param = XYparameter[vec['from']] 
            B = param['B']
            Km = param['Km']
            n =  param['n']
            if vec['effect'] == -1:
                addTerm += repressor(B,Km,n,inputVal)
            elif vec['effect'] == 1:
                addTerm += activator(B,Km,n,inputVal)
            else:
                addTerm += 0
    answer = firstTerm+addTerm
    return answer

# B: production rate of protein
# K: repression coefficient, higher value requires more repressor to repress expression
# n: hill coefficient, cooperativity
#   - higher n causes system to behave as if input is digital (1/0)
# inputVal: concentration of repressor protein 
def repressor(B,K,n,inputVal):
    '''
        Repressor hill equation
    '''
    try:
        return B/(1+(inputVal/K)**n)
    # If K is zero, should never occur (K > 0)
    except ZeroDivisionError:
        try:
            return B*(K**n)/(inputVal**n+K**n)
        # If both K and Input are zero
        except ZeroDivisionError:
            return 0

# B: production rate of protein
# K: activation coefficient, higher value requires more activator to activate expression
# n: hill coefficient, higher value means more cooperativity
#   - higher n causes system to behave as if input is digital (1/0)
# inputVal: concentration of activator protein 
def activator(B,K,n,inputVal):
    '''
        Activator hill equation
    '''
    try:
        return B/(1+(K/inputVal)**n)
    # If Input is zero
    except ZeroDivisionError:
        try:
            return B*(inputVal**n)/(inputVal**n+K**n)
        # If both K and Input are zero
        except ZeroDivisionError:
            return 0
            
# File read should contain json format of initial gene conditions, parameters, and vectors
fileLoc = 'masterXYInputFile1.json'

# Code for plotting configurations of X/Y nodes interactions
xMax = 100
itr = 10000
t = np.linspace(0, xMax, itr)
X_Y_initCond, X_Y_params, X_Y_vectors = getFromJsonFile(fileLoc)
firstState = []
keys = X_Y_initCond.keys()
keys.sort()
# The initial conditions should correspond to the sorted keys
for key in keys:
    firstState.append(X_Y_initCond[key])
    
for i in range(len(X_Y_vectors)):
    # firstState maps to state when function f is called in the odeint function
    # index i for X_Y_vectors gets us an individual array of objects
    soln_f = odeint(f, firstState, t, args = (X_Y_params, X_Y_vectors[i]))
    X = soln_f[:,0]
    Y = soln_f[:,1]
    
    # Titles graph according to vectors for clarity    
    title = ''
    for j in range(4):
        title += 'from ' + str(X_Y_vectors[i][j]['from']) +' to '+ str(X_Y_vectors[i][j]['to']) + ": " + str(X_Y_vectors[i][j]['effect']) + ' ' 
    title = title[:-2]

    # Plots the transcriptional dynamics
    plt.figure()
    #plt.axis([0, 100, 0 , 9000])
    plt.plot(t,X, 'r-', label = 'X')
    plt.plot(t,Y, 'b-', label = 'Y')
    plt.xlabel('Time (min)')
    plt.ylabel('Proteins per cell')
    plt.title(title)
    plt.legend(loc=0)
