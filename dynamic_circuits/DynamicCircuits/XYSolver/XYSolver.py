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

def f(state, t, parameter, vector):
    #Set minimum value to 0.0 if negative.
    for i in range(len(state)):
        if state[i]<0 :
            state[i] = 0.0
            
    #Concentration of proteins at this state
    initProtein = state
    #Initialize the lists that hold the change in values (deltaValues).
    dinitProtein_dt = []

    #For each gate we want to calculate the change in its mRNA and protein
    for i in range(len(state)):
        #Retrieve some current Properties
        initProteinVal = initProtein[i]
        #Get the changes and add it to the list holding the changes
        if i == 0:
            inp = 'X'
        else:
            inp = 'Y'
        dinitProtein_dt.append(getProteinChange(initProtein,parameter[inp],vector))
    return dinitProtein_dt
    
def getFromJsonFile(fileLoc):
    '''
    Takes in a file location and returns the values necessary to generate a graph
    '''
    myFile = open(fileLoc,'r')
    data = json.load(myFile)
    X_Y_params = data[0]
    X_Y_vectors = data[1:]
    return X_Y_params, X_Y_vectors

# initProtein is a dictionary, parameter is a dictionary of 2 dictionaries, 
# vector is a single list of dictionaries of 'to's & 'from's
def getProteinChange(initProtein, param, vector):
    # Degradation term and basal rate
    firstTerm = -param['a']*initProtein + param['B0']
    name = param['name']
    B = param['B']
    Km = param['Km']
    n = param['n']

    # Regulation terms
    addTerm = 0
    
    # Activity of node X
    if name == 'X':
        for vec in vector:
            if vec['to'] == 'X':
                if vec['from'] == 'X':
                    inputVal = initProtein
                elif vec['from'] == 'Y':
                    inputVal = initProtein
                if vec['effect'] == -1:
                    addTerm += repressor(B,Km,n,inputVal)
                elif vec['effect'] == 1:
                    addTerm += activator(B,Km,n,inputVal)
                else:
                    addTerm += 0
        answer = firstTerm+addTerm
    elif name == 'Y':
        for vec in vector:
            if vec['to'] == 'Y':
                if vec['from'] == 'X':
                    inputVal = initProtein['X']
                elif vec['from'] == 'Y':
                    inputVal = initProtein['Y']
                if vec['effect'] == -1:
                    addTerm += repressor(B,Km,n,inputVal)
                elif vec['effect'] == 1:
                    addTerm += activator(B,Km,n,inputVal)
                else:
                    addTerm += 0
        answer = firstTerm+addTerm
    return answer

#Repressor and Activator Hill Equations
def repressor(B,Km,n,Input):
    '''
        Repressor hill equation
    '''
    try:
        return B/(1+(Input/Km)**n)
    #if Km is zero
    except ZeroDivisionError:
        try:
            return B*(Km**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0

def activator(B,Km,n,Input):
    '''
        Activator hill equation
    '''
    try:
        return B/(1+(Km/Input)**n)
    #if Input is zero
    except ZeroDivisionError:
        try:
            return B*(Input**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0


