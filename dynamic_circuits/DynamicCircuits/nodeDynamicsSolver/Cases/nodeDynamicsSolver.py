# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
Program that takes in a json file of an array containing 2 dictionaries of the 
initial state of the system of nodes and gene parameters as well as an array 
of objects specifying the regulation motifs. (See masterNodeInputFile.json)

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
    Reads in a json file directory and returns values needed for plotting.
    '''
    myFile = open(fileLoc,'r')
    data = json.load(myFile)
    node_inic = data[0] # a dictionary of initial conditions
    node_params = data[1] # a dictionary containing 2 dictionaries of protein parameters
    node_vectors = data[2:] # a list of arrays containing objects defining regulation motifs
    return node_inic, node_params, node_vectors

# state: the initial condition at which the differential would be determined
# t: the list of time values
# parameter: a dictionary of gene parameters where the keys are the gene names
#   -For example, 'X' would be a key to a dictionary of parameters (see masterNodeInputFile.json)
# vectors: an array of objects with each specifying a regulation motif  
def f(state, t, parameter, vectors):
    '''
    Function passed through odeint() that represents the system of differential 
    equations of the rate of change of each protein.
    '''
    proteinNames = parameter.keys()
    proteinNames.sort()
    initProtein = {}
    
    # By convention, the first index of state corresponds to the first key in the initProtein
    # Sets the minimum value to 0.0 if the concentration of protein at this state is negative.
    for i in range(len(state)):
        if state[i]<0 :
            state[i] = 0.0
        initProtein[proteinNames[i]] = state[i]
            
    # Initializes the lists that hold the change in values (deltaValues).
    dinitProtein_dt = []

    # For each protein we want to calculate the change in its levels
    for name in proteinNames:
        # Gets the changes and adds it to the list holding the changes
        dinitProtein_dt.append(getProteinChange(name,initProtein,parameter,vectors))
    return dinitProtein_dt

# proteinName: a string of the variable that rate of change is being measured
# initProtein: a dictionary of the initial conditions, 
# parameter: a dictionary of 2 dictionaries of protein parameters, 
# vectors: a single array of objects representing regulation motifs.
def getProteinChange(proteinName, initProtein, parameters, vectors):
    param = parameters[proteinName]

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
            param = parameters[vec['from']] 
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
            

# json file should contain an array of 2 objects and an array representing
# initial gene conditions, parameters, and vectors, respectively.
jsonLoc = 'repressilator.json'
# png file name & location where png will be saved
pngLoc = 'repressilator.png'

# Code for plotting configurations of X/Y nodes interactions
xMax = 1000
itr = 5*xMax
t = np.linspace(0, xMax, itr)
node_initCond, node_params, node_vectors = getFromJsonFile(jsonLoc)
firstState = []
proteinNames = node_initCond.keys()
numProteins = len(proteinNames)
proteinNames.sort()
# The initial conditions should correspond to the sorted keys
for key in proteinNames:
    firstState.append(node_initCond[key])
    
for i in range(len(node_vectors)):
    # firstState maps to state when function f is called in the odeint function
    # index i for node_vectors gets us an individual array of objects
    soln_f = odeint(f, firstState, t, args = (node_params, node_vectors[i]))

    # Plots the transcriptional dynamics
    plt.figure()
    for j in range(numProteins): 
        plt.plot(t,soln_f[:,j], label = proteinNames[j])
    plt.xlabel('Time (min)')
    plt.ylabel('Proteins per cell')
    plt.title('Repressilator')
    plt.legend(loc=0)
    #Saves file as png file
    plt.savefig(pngLoc)
    
    
    
    
    
    
    
    