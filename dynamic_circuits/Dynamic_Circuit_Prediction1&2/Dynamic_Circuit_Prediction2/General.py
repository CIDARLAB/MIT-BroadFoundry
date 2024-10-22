# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 22:41:15 2015

Does the same as Dynamic_Circuit_Prediction but uses REU values and Pmin and Pmax
values and ignores mRNA, only focusing on protein actions.

@author: Arinze
"""

import matplotlib.pyplot as plt
import update
import GeneralJsonIO
from scipy.integrate import odeint
import numpy as np
import DifferentialSolver
import time

#input is a json file
def generateDynamicCircuitGraphs(fileName, makeBarGraphs, makeOtherGraphs, isSequential): 
    startTime = time.time()
#---------------------------------------------------------------
    #Get all the values from the given file.
    #input_and_logic_gate_dictionaries: dictionary for each gate with all properties
    #properties: 
    time_axis_params,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names = GeneralJsonIO.getFromJsonFile(fileName)
#---------------------------------------------------------------

    #Generate the time grid
    numInputs = len(input_and_logic_gate_names)-len(logic_gate_names)
    
    #for example: 0, 500, 10,000: [0,....500] (9,998 points in between)
    t = np.linspace(0, time_axis_params['xMax'], time_axis_params['itr'])
    numGates = len(logic_gate_names)

    initc = []
    for i in range(numGates):
        #Add the initial REU in the order they appear in logic_gate_names.
        initc.append(input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]['REUi'])

    #Use odeint to calculate the concentrations over the time steps

    #arguments:
    #f: function.  takes state and a time point, returns the delta of each value.  First 'state' is equal to initc.
    #... adds deltas to state to calculate the next state values.
    #initc: the true initial condition for the very beginning.
    #t: coarse, odeint breaks it into finer timepoints.
    #args: every other input that f takes that is not state or t
    soln_f = odeint(f, initc, t, args=(input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,numGates),hmax=350)
    #soln_f = DifferentialSolver.differentialSolver(f, initc, t, args=(input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,numGates))
    
    #soln_f is a list of length len(t) containing lists of length len(initc).
    #The ith index of soln_f corresponds to the protein concentration at the
    #time of the ith index of t.
    
    REU_f = []
    for i in range(numGates):
        REU_f.append(soln_f[:,i])
    #REU_f is a list of length len(initc) where each list is of length len(t)
    #and each list contains the concentrations of the proteins for a gate at
    #a given time
        
    #These are the possible graph lines/markers. This is used to make sure each
    #line on a graph is unique.
    color = ['b','g','r','c','m','y','k']
    linestyle = ['--','-','-.',':']
    marker = ['','.',',','o','v','^','<','>','s','p','*','h','+','x','D','|','_']
    
    if makeOtherGraphs:
        #Visualizing Inputs
        #Make the graphs for the inputs
        inputs_f = []
        inputNames = []
        for thing in input_and_logic_gate_names:
            #For each element in input_and_logic_gate_names, if it is not in logic_gate_names, it is an input.
            if thing not in logic_gate_names:
                #Save the name and a list of its values for each time step.
                inputNames.append(thing)
                
                #see also GeneralJSONIO.py... gets the input function and params
                #func[0] is the input function
                func = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(thing)]['INPUT']
                temp = []
                
                #calculate what the input value is for each timepoint
                for timepos in t:
                    temp.append(func[0](timepos,*func[1:]))
                inputs_f.append(temp)

        #set plot ranges
        ymax=1
        for i in inputs_f:
            m = max(i)*10
            if m>ymax:
                ymax = m
                
        #make one plot that has all the input values (may be overlapping)
        plt.figure()
        plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymax=ymax, ymin=0.1)  
        #Plot each input against time with its name as the label
        for i in range(len(inputNames)):
            tempInput = list(inputs_f[i])
            #tempInput.reverse()
            plt.plot(t,tempInput,color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label=inputNames[i])
        plt.xlabel('Time (min)')
        plt.yscale('log')
        plt.ylabel('Input molecules per cell')
        plt.title('Concentration of Inputs Over Time')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
            
        #Visualizing REU (each gate on one graph, except inputs/outputs)
        plt.figure()
        #Plot each REU against time with its name as the label
        ymax = 1
        for i in range(numGates):
            tempREU = list(REU_f[i])
            tempMax = max(tempREU)
            if tempMax>ymax:
                ymax = tempMax
            plt.plot(t,tempREU,color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label = logic_gate_names[i])
        plt.xlabel('Time (min)')
        plt.yscale('log')
        plt.ylabel('REU')
        plt.title('REU Over Time')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ymax = ymax *10
        plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymin=0.1, ymax = ymax)
        plt.show()
    
        #Makes a graph that only contains the outputs.
        print "Everything but outputs stripped from graph."
        #Visualizing REU
        plt.figure()
        #Plot each output REU against time with its name as the label
        ymax = 1
        for i in range(numGates):
            tempREU = list(REU_f[i])
            tempMax = max(tempREU)
            if tempMax>ymax:
                ymax = tempMax
            tempGateType = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]["TYPE"]
            if tempGateType=="OR" or tempGateType=="BUFFER":
                plt.plot(t,tempREU,color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label = logic_gate_names[i])
        plt.xlabel('Time (min)')
        plt.yscale('log')
        plt.ylabel('Transcripts per cell')
        plt.title('REU Over Time')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ymax = ymax *10
        plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymin=0.1, ymax = ymax)
        plt.show()

#-------------------------------------------------------------------------------
    if isSequential:
        return # don't bother making bar-graphs or getting scores
    """    
    Make bar graphs for each genes mRNA and protein. It will show its mRNA(or protein)
    at each truth value. Bars that are expected to be high are blue. Bars that
    are expected to be low are red.
    """
    REUListOfListOfValues = []
    #Tries to estimate the value that the protein converges or mRNA 
    #converges at for each gate.
    for i in range(numGates):
        REUListOfListOfValues.append(getListOfConvergingValues(numInputs, time_axis_params['itr'], REU_f[i]))
    #Get the binary values for the bar graph x axis labels
    numBinaries = 2**numInputs
    binVals = getBinaryInOrder(numBinaries)
    t = np.arange(numBinaries)
    #Constant that ensures the x axis labels are properly spaced
    width = 0.425

    scoreDict = {}
    
    #Make a bar graph for each gene
    for i in range(numGates):
        REUVals = list(REUListOfListOfValues[i])
        REUymax = max(REUVals)*10
        currGateProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]
        tv = currGateProperties['EXPECTED_PROTEIN']
        name = str(currGateProperties['NAME'])
        if makeBarGraphs:
            #Make the bar plot for the REU
            plt.figure()
            plt.axis(xmin=0, xmax=numBinaries, ymin=0.1, ymax=REUymax)
            barlist = plt.bar(t,REUVals)
            #If the expected is low, change the color of the bar to red.
            for j in range(numBinaries):
                if tv[j]=="0":
                    barlist[j].set_color('r')
            #Space the labels properly
            plt.xticks(t+width, binVals)  
            plt.xlabel('TruthValues')
            plt.yscale('log')
            plt.ylabel('REU')
            plt.title('REU for ' + name)
            plt.show()
        lowestHigh = REUymax
        highestLow = 0.0
        for i in range(numBinaries):
            tempVal = REUVals[i]
            if tv[i] == "1" and tempVal<lowestHigh:
                lowestHigh = tempVal
            elif tv[i] == "0" and tempVal>highestLow:
                highestLow = tempVal
        
        gateScore = lowestHigh/highestLow
        scoreDict[name] = gateScore
        
#    print "Gate scores:"
#    print scoreDict
    endTime = time.time()
#    print round(endTime-startTime,2), "seconds"
    return scoreDict
    
def getBinaryInOrder(maxVal):
    """
    returns a list of binary numbers (as strings) from 0 to maxVal with
    the length of each binary number equal to the length of maxVal as a binary.
    this is done by adding leading zeros to the binary
    """
    binVals = []
    #get a list of the numbers in decimal
    decVals = range(maxVal)
    #reverse the list so the largest number is first
    decVals.reverse()
    for i in decVals:
        #convert to binary
        temp="{0:b}".format(i)
        #if it is the first item then add it to the list it will be the size 
        #that all of the binary values should be because the first item will be
        #(11...1) where as other ones will be missing the leading zeros.
        if len(binVals) == 0:
            binVals.append(temp)
        
        else:
            #Add zeroes to the front until the value is the correct size
            while len(temp)<len(binVals[0]):
                temp = "0"+temp
            binVals.append(temp)
    #Reverse the order so the smallest value is back in front
    binVals.reverse()
    return binVals

def f(state, t, input_and_logic_gate_dictionaries, input_and_logic_gate_names, logic_gate_names, numGates):
    """
    takes in a current state and a time point (along with other arguments)
    and returns the slope of the time vs state graph at the given time point.
    """
    for i in range(len(state)):
        if state[i]<0 :
            state[i] = 0.0
            
    #results.length equals the number of gates (including outputs, excluding inputs)
    #values in results are dX/dt, dY/dt, etc...
    results = []
    
    #For each gate we want to calculate the change in its REU
    for i in range(numGates):
        
        #Retrieve some current Properties
        gateProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]
        initREUVal = state[i]
        #Get the changes and add it to the list holding the changes
        results.append(update.getREUChange(t,initREUVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,state))
        
    #return an array of slopes
    #odeint will multiply the slope by the time step to get the delta
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
    return output
'''
_______________________________________________________________________________
'''

            