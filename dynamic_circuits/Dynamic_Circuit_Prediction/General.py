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
def generateDynamicCircuitGraphs(fileName, makeBarGraphs): 
#---------------------------------------------------------------
    #Get all the values from the given file.
    print "Loading values from file."
    time_axis_params,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names = GeneralJsonIO.getFromJsonFile(fileName)
    print "Finished loading values from file."
#---------------------------------------------------------------

    #Generate the time grid
    numInputs = len(input_and_logic_gate_names)-len(logic_gate_names)
    t = np.linspace(0, time_axis_params['xMax'], time_axis_params['itr'])
    #This generates a list of every time point at which the input 
    #dramatically changes. We want odeint to be more careful at these points
    carefult = np.linspace(0, time_axis_params['xMax'], numInputs+2) 
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
    soln_f = odeint(f, initc, t, args=(input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names),tcrit=carefult)
    
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
    ymax=1
    for i in inputs_f:
        m = max(i)*10
        if m>ymax:
            ymax = m
    plt.figure()
    plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymax=ymax, ymin=1)  
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
    
    #Visualizing mRNA concentrations
    plt.figure()
    #Plot each mRNA against time with its name as the label
    ymax = 1
    for i in range(numGates):
        tempmRNA = list(mRNA_f[i])
        #tempmRNA.reverse()
        tempMax = max(tempmRNA)
        if tempMax>ymax:
            ymax = tempMax
        plt.plot(t,tempmRNA,color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label = logic_gate_names[i])
    plt.xlabel('Time (min)')
    plt.yscale('log')
    plt.ylabel('Transcripts per cell')
    plt.title('Concentration of mRNA Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ymax = ymax *10
    plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymin=1, ymax = ymax)
    plt.show()
    
    #Visualizing protein concentrations
    plt.figure()
    #Plot each protein against time with its name as the label
    ymax = 1
    for i in range(numGates):
        tempProtein = list(protein_f[i])
        #tempProtein.reverse()
        tempMax = max(tempProtein)
        if tempMax>ymax:
            ymax = tempMax
        plt.plot(t,tempProtein,color=color[i%7],linestyle=linestyle[(i/7)%4],marker=marker[((i/7)/4)%17],label = logic_gate_names[i])
    plt.xlabel('Time (min)')
    plt.yscale('log')
    plt.ylabel('Proteins per cell')
    plt.title('Concentration of Protein Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ymax = ymax *10
    plt.axis(xmin=time_axis_params['xMin'], xmax=time_axis_params['xMax'], ymin=1, ymax = ymax)
    plt.show()
    
#-------------------------------------------------------------------------------
    if not makeBarGraphs:
        return
    """    
    Make bar graphs for each genes mRNA and protein. It will show its mRNA(or protein)
    at each truth value. Bars that are expected to be high are blue. Bars that
    are expected to be low are red.
    """
    proteinListOfListOfValues = []
    mRNAListOfListOfValues = []
    #Tries to estimate the value that the protein converges or mRNA 
    #converges at for each gate.
    for i in range(numGates):
        proteinListOfListOfValues.append(getListOfConvergingValues(numInputs, time_axis_params['itr'], protein_f[i]))
        mRNAListOfListOfValues.append(getListOfConvergingValues(numInputs, time_axis_params['itr'], mRNA_f[i]))    
    #Get the binary values for the bar graph x axis labels
    numBinaries = 2**numInputs
    binVals = getBinaryInOrder(numBinaries)
    t = np.arange(numBinaries)
    #Constant that ensures the x axis labels are properly spaced
    width = 0.425

    scoreDict = {}
    
    #Make two bar graphs for each gene. One for mRNA and one for protein
    for i in range(numGates):

        proteinVals = list(proteinListOfListOfValues[i])
        mRNAVals = list(mRNAListOfListOfValues[i])

        #Add 10% so the ymax is slightly higher than the highest bar
        mrnaymax = max(mRNAVals)*10
        proteinymax = max(proteinVals)*10
        #Get the properties of the gate so we can find the name and 
        #expected highs and lows
        currGateProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(logic_gate_names[i])]
        tv = currGateProperties['EXPECTED_PROTEIN']
        name = str(currGateProperties['NAME'])
        
        #Make the bar plot for the mRNA
        plt.figure()
        plt.axis(xmin=0, xmax=numBinaries, ymin=1, ymax=mrnaymax)
        barlist = plt.bar(t,mRNAVals)
        #If the expected is low, change the color of the bar to red.
        for j in range(numBinaries):
            if tv[j]=="0":
                barlist[j].set_color('r')
        #Space the labels properly
        plt.xticks(t+width, binVals)  
        plt.xlabel('TruthValues')
        plt.yscale('log')
        plt.ylabel('mRNA per cell')
        plt.title('mRNA Levels for ' + name)
        plt.show()
        
        #Make the bar plot for the protein
        plt.figure()
        plt.axis(xmin=0, xmax=numBinaries, ymin=1, ymax=proteinymax)
        barlist = plt.bar(t,proteinVals)
        #If the expected is low, change the color of the bar to red.
        for j in range(numBinaries):
            if tv[j]=="0":
                barlist[j].set_color('r')
        #Space the labels properly
        plt.xticks(t+width, binVals)  
        plt.xlabel('TruthValues')
        plt.yscale('log')
        plt.ylabel('Proteins per cell')
        plt.title('Protein Levels for ' + name)
        plt.show()
        
        #Calculate the gate score
        lowestHigh = proteinymax
        highestLow = 0.0
        for i in range(numBinaries):
            tempVal = proteinVals[i]
            if tv[i] == "1" and tempVal<lowestHigh:
                lowestHigh = tempVal
            elif tv[i] == "0" and tempVal>highestLow:
                highestLow = tempVal
        gateScore = lowestHigh/highestLow
        scoreDict[name] = gateScore
        
    print "Gate scores:"
    print scoreDict
    
def getBinaryInOrder(maxVal):
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
    #print output
    return output
'''
_______________________________________________________________________________
'''
