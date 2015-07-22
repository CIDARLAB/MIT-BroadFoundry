# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 22:41:15 2015

@author: Arinze
"""

import matplotlib.pyplot as plt
import inputs
import update
import GeneralJsonIO
from scipy.integrate import odeint
import numpy as np

#input is a json file
def generateDynamicCircuitGraphs(fileName): 
#---------------------------------------------------------------
    #Get all the values from the given file.
    print "Loading values from file."
    time_axis_params,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names = GeneralJsonIO.getFromJsonFile(fileName)
    print "Finished loading values from file."
#---------------------------------------------------------------
    #Generate the time grid
    t = np.linspace(0, time_axis_params['xMax'], time_axis_params['itr']) 
    #Initialize the list of mRNA initial concentrations and protein initial concentrations.
    mRNAi = []
    proteini = []
    for i in range(len(logic_gate_names)):
        #If a degredation constant alpha was not given, calculate it from the half-life
        if input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['am']=='NONE':
            input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['am'] = np.log(2)/input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['m_HALFLIFE']
        if input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['ap']=='NONE':
            input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['ap'] = np.log(2)/input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['p_HALFLIFE']
        #Add the initial mRNA and protein concentrations in the order they appear
        #in logic_gate_names.
        mRNAi.append(input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['Mi'])
        proteini.append(input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['Pi'])
    #Combine into one vector of starting concentrations
    initc = mRNAi + proteini
    #Use odeint to calculate the concentrations over the time steps
    soln_f = odeint(f, initc, t, args=(input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names))
    
    #Separate the mRNA and protein concentration vectors. The order they appear
    #in the list should be the same order that they are in in logic_gate_names.
    mRNA_f = []
    protein_f = []
    numGates = len(logic_gate_names)
    for i in range(numGates):
        mRNA_f.append(soln_f[:,i])
        protein_f.append(soln_f[:,numGates+i])
        
    #These are the possible graph lines/markers.
    color = ['b','g','r','c','m','y','k']
    linestyle = ['--','-','-.',':']
    marker = ['','.',',','o','v','^','<','>','s','p','*','h','+','x','D','|','_']
      
    #Visualizing Inputs
    inputs_f = []
    inputNames = []

#-------------------------------------------------------------------------------
    numInputs = len(input_and_logic_gate_names)-len(logic_gate_names)
    proteinListOfListOfValues = []
    mRNAListOfListOfValues = []
    for i in range(numGates):
        proteinListOfListOfValues.append(getListOfConvergingValues(numInputs, time_axis_params['itr'], protein_f[i]))
        mRNAListOfListOfValues.append(getListOfConvergingValues(numInputs, time_axis_params['itr'], mRNA_f[i]))    
    
    print proteinListOfListOfValues[0]
    print mRNAListOfListOfValues[0]
#-------------------------------------------------------------------------------
   
    for thing in input_and_logic_gate_names:
        #For each element in input_and_logic_gate_names, if it is not in logic_gate_names, it is an input.
        if thing not in logic_gate_names:
            #Save the name and a list of its values for each time step.
            inputNames.append(thing)
            func = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(thing)]['INPUT']
            temp = []
            for timepos in t:
                temp.append(func[0](timepos,*func[1:]))
            inputs_f.append(temp)
    ymax=0
    for i in inputs_f:
        m = max(i)*1.1
        if m>ymax:
            ymax = m
    plt.figure()
    plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymax=ymax, ymin=-10)
   
   #Plot each input against time with its name as the label
    for i in range(len(inputNames)):
        plt.plot(t,inputs_f[i],color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label=inputNames[i])
    plt.xlabel('Time (min)')
    plt.ylabel('Input molecules per cell')
    plt.title('Concentration of Inputs Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
    #Visualizing mRNA concentrations
    plt.figure()
    plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'])
    #Plot each mRNA against time with its name as the label
    for i in range(numGates):
        plt.plot(t,mRNA_f[i],color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label = logic_gate_names[i])
    plt.xlabel('Time (min)')
    plt.ylabel('Transcripts per cell')
    plt.title('Concentration of mRNA Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    #Visualizing protein concentrations
    plt.figure()
    plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'])
    #Plot each protein against time with its name as the label
    for i in range(numGates):
        plt.plot(t,protein_f[i],color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label = logic_gate_names[i])
    plt.xlabel('Time (min)')
    plt.ylabel('Proteins per cell')
    plt.title('Concentration of Protein Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

def f(state, t, input_and_logic_gate_dictionaries, input_and_logic_gate_names, logic_gate_names):
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
'''
Alex was working here. Code for creating a list of the values that the graph 
will converge to given for mRNA and proteins
_______________________________________________________________________________
'''
    #List of truthvalues of gates
    #tvs = input_and_logic_gate_dictionaries['TRUTHVALUES']    

#Function for obtaining list of values given a graph of the output graph
def getListOfConvergingValues(numInputs, numItr, opGraph):
    output = []    
    copy = opGraph
    interval = numItr/2**numInputs
    for i in range(2**numInputs):
        output.append(copy[(i+1)*interval-2])
    print output
    return output
'''
_______________________________________________________________________________
'''
