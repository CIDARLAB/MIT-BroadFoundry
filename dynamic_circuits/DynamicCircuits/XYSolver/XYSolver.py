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
'''
def f(state, t):
        X = state[0]
        Y = state[1]
        
        for i in state:
        # Repressor Protein Equations
        
        if(s == 'null'):
            dX_dt = B0X + BX - aX*X # Nothing
            dX_dt = B0X + BX/(1 + np.power(X/Km,n)) - aX*X # Auto-Repression
            dX_dt = B0X + BX*np.power(X,n)/(np.power(X,n)+np.power(Km,n)) - aX*X # Auto-Activation
            
            dX_dt = B0X + BX - aX*X # Nothing
            dX_dt = B0X + BX * np.power(Y,n)/(np.power(Y,n)+np.power(Km,n)) - aX*X # Activation
            dX_dt = B0X + BX/(1 + np.power(Y/Km,n)) - aX*X # Repression
            
            
            dY_dt = B0Y + BY - aY*Y
            dY_dt = B0Y + BY*np.power(X,n)/(np.power(X,n)+np.power(Km,n)) - aY*Y
            dY_dt = B0Y + BY/(1 + np.power(X/Km,n)) - aY*Y
        elif():
        return [dX_dt, dY_dt]
'''
def f(state, t, input_and_logic_gate_dictionaries, input_and_logic_gate_names, logic_gate_names):
    for i in range(len(state)):
        if state[i]<0 :
            state[i] = 0.0
            
    #Concentration of mRNA and protein at this state
    numGates = len(logic_gate_names)
    #Split mRNA values from protein values
    initmRNA = state[:numGates]
    initProtein = state[numGates:]
    #Initialize the lists that hold the change in values (deltaValues).
    dinitmRNA_dt = []
    dinitProtein_dt = []

    #For each gate we want to calculate the change in its mRNA and protein
    for i in range(numGates):
        #Retrieve some current Properties
        gateProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]
        initmRNAVal = initmRNA[i]
        initProteinVal = initProtein[i]
        #Get the changes and add it to the list holding the changes
        dinitmRNA_dt.append(update.getmRNAChange(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein))
        dinitProtein_dt.append(update.getProteinChange(initProteinVal,initmRNAVal,gateProperties))
    #recombine into one list
    results = dinitmRNA_dt + dinitProtein_dt
    return results
    
def getFromJsonFile(fileLoc):
    '''
    Takes in a file location and returns the values necessary to generate a graph
    '''
    myFile = open(fileLoc,'r')
    data = json.load(myFile)
    X_Y_params = data[0]
    X_Y_vectors = data[1:]
    return X_Y_params, X_Y_vectors

def getProteinChange(t, initProtein, parameter, vector):
    if something == 'X':
        param = parameter['X']
        param2 = parameter['Y']
    elif something == 'Y':
        param = parameter['Y']
        param2 = parameter['X']
    # Degradation term and basal rate
    firstTerm = -param['a']*initProtein + param['B0']

    B = param['B']
    Km = param['Km']
    n = param['n']

    # Autoregulation term
    secondTerm = 0
    # Regulation term
    thirdTerm = 0

    
    #Repression
    if vector['autoX'] == -1:
        Input = initProtein['X']    #might not be instantiated yet.
        secondTerm += repressor(B,Km,n,Input)
    elif vector['autoX'] == 0:
        secondTerm += 0
    elif vector['autoX'] == 1:
        secondTerm += activator(B,Km,n,Input)
time_axis_params,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names = GeneralJsonIO.getFromJsonFile(fileName)


def getmRNAChange(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    #first term is degradation rate other terms are from
    firstTerm = -gateProperties['am']*initmRNAVal
    
    #These gates only have one input so they only have two terms
    if gateProperties['TYPE']=='BUFFER' or gateProperties['TYPE']=='NOT':
        #Second term
        mB = gateProperties['mB']
        Km = gateProperties['Km']
        n = gateProperties['n']
        InputName = gateProperties['INPUT']
        InputProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName)]
        #If this is another gate then just get the concentration of the protein from it.
        if InputProperties['TYPE'] != 'INPUT':
            Input = initProtein[logic_gate_names.index(InputName)]
        #If it is an input, then use the equation specified for the input to determine
        #the concentration of the input at that time step.
        else:
            Inputfunc = InputProperties['INPUT']
            Input = Inputfunc[0](t,*Inputfunc[1:])
        if gateProperties['INPUT_EFFECT']=="REPRESS":
            secondTerm = repressor(mB,Km,n,Input)
        elif gateProperties['INPUT_EFFECT']=="ACTIVATE":
            secondTerm = activator(mB,Km,n,Input)
        answer =  firstTerm + secondTerm   
        
    #These gates have two inputs so they will have a total of three terms
    elif gateProperties['TYPE']=='OR' or gateProperties['TYPE']=='NOR':
        #Second term
        mB1 = gateProperties['mB1']
        Km1 = gateProperties['Km1']
        n1 = gateProperties['n1']
        InputName1 = gateProperties['INPUT1']
        InputProperties1 = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName1)]
        #If this is another gate then just get the concentration of the protein from it.
        if InputProperties1['TYPE'] != 'INPUT':
            Input1 = initProtein[logic_gate_names.index(InputName1)]
        #If it is an input, then use the equation specified for the input to determine
        #the concentration of the input at that time step.
        else:
            Inputfunc1 = InputProperties1['INPUT']
            Input1 = Inputfunc1[0](t,*Inputfunc1[1:])
        if gateProperties['INPUT1_EFFECT']=="REPRESS":
            secondTerm = repressor(mB1,Km1,n1,Input1)
        elif gateProperties['INPUT1_EFFECT']=="ACTIVATE":
            secondTerm = activator(mB1,Km1,n1,Input1)
        #Third Term
        mB2 = gateProperties['mB2']
        Km2 = gateProperties['Km2']
        n2 = gateProperties['n2']
        InputName2 = gateProperties['INPUT2']
        InputProperties2 = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName2)]
        #If this is another gate then just get the concentration of the protein from it.
        if InputProperties2['TYPE'] != 'INPUT':
            Input2 = initProtein[logic_gate_names.index(InputName2)]
        #If it is an input, then use the equation specified for the input to determine
        #the concentration of the input at that time step.
        else:
            Inputfunc2 = InputProperties2['INPUT']
            Input2 = Inputfunc2[0](t,*Inputfunc2[1:])
        if gateProperties['INPUT2_EFFECT']=="REPRESS":
            thirdTerm = repressor(mB2,Km2,n2,Input2)
        elif gateProperties['INPUT2_EFFECT']=="ACTIVATE":
            thirdTerm = activator(mB2,Km2,n2,Input2)
        answer = firstTerm + secondTerm + thirdTerm        
        
    return answer
    

#Repressor and Activator Hill Equations
def repressor(mB,Km,n,Input):
    '''
        Repressor equation
    '''
    try:
        return mB/(1+(Input/Km)**n)
    #if Km is zero
    except ZeroDivisionError:
        try:
            return mB*(Km**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0
def activator(mB,Km,n,Input):
    '''
        Activator equation
    '''
    try:
        return mB/(1+(Km/Input)**n)
    #if Input is zero
    except ZeroDivisionError:
        try:
            return mB*(Input**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0