# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:17:06 2015

@author: Arinze
"""
import json
import Gate
import itertools
import OptimalCircuit
import General
import copy
import random
import time
import matplotlib.pyplot as plt

inputsDir = "JsonFiles/Libraries/InputLibrary.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary.json"
Libraries = [inputsDir, repressorsDir, outputsDir]
truthValueExampleFileLoc = "JsonFiles/SplitByTruthValue_OR/01101001.json" #01101001
truthValueExampleFileLoc2 = "JsonFiles/SplitByTruthValue_OR/11110010.json" #11110010
placeToSaveTruthValueExample = "JsonFiles/01101001Graph.json"
exampleNetlist = ['NOR(W1,a,b)','NOR(W2,W1,c)','NOR(W3,a,c,d)','NOR(W4,W1,a,d)','NOR(W5,W2,W3,W4)']

def findOptimalAssortmentHill(netlist,Libraries=Libraries, smallestScoreAllowed=10,numTraj=5):
    #Keep track of how long it takes to find each improvement.
    xvals = [0]
    yvals = [0]
    startTime = time.time()
    #random.seed(0)
    
    #this is where the graph will be saved as a json file
    scratchFile = "JsonFiles/optimizingScratchWorkHill.json"
    #make the starting graph with the default values
    graph, allGates, Inputs, Repressors, Outputs = OptimalCircuit.makeGraphFromNetlist(netlist, True)
    #if it is a sequential circuit this process will not work
    if graph.hasLoop():
        print "This circuit is sequential and cannot be scored."
        return
    
    #load in all the gates from the libraries
    allInputs = makeInputs(Libraries[0])
    allRepressors = makeRepressors(Libraries[1])
    allOutputs = makeOutputs(Libraries[2])
    numInputs = len(Inputs)
    numRepressors = len(Repressors)
    numOutputs = len(Outputs)
    chosenInputs = selectInputs(allInputs,numInputs)
    chosenOutputs = selectOutputs(allOutputs,numOutputs)
    for i in range(numInputs):
        graph.swapGates(chosenInputs[i],Inputs[i])
    outputNames = []
    for i in range(numOutputs):
        graph.swapGates(chosenOutputs[i], Outputs[i])
        outputNames.append(chosenOutputs[i].getName())

    #initialize the global best
    globalBestGraph = None
    globalBestScore = 0
    globalBestScoreDict = None
    #get all possible swaps
    standardAllRepressorSwappingPairs = list(itertools.combinations(allRepressors,2))
    
    #start the trajectory
    for i in range(numTraj):
        #randomly assign the starting point
        chosenStartRepressors = selectStartingRepressors(allRepressors,numRepressors)
        for j in range(numRepressors):
            graph.swapGates(chosenStartRepressors[j], Repressors[j])
        Repressors = chosenStartRepressors
        while graph.hasRepeatedRepressor():
            chosenStartRepressors = selectStartingRepressors(allRepressors,numRepressors)
            for j in range(numRepressors):
                graph.swapGates(chosenStartRepressors[j], Repressors[j])
            Repressors = chosenStartRepressors

        #initialize the local best
        localBestGraph = None
        localBestScore = 0
        localBestScoreDict = None
        #get the score of this starting point
        graph.writeToJson(scratchFile)
        firstScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
        firstResults = getScoreFromDict(outputNames, firstScoreDict, smallestScoreAllowed)
        currScore = firstResults[0]
        if firstResults[1]:
            localBestGraph = copy.deepcopy(graph)
            localBestScore = firstResults[0]
            localBestScoreDict = firstScoreDict
            if localBestScore>globalBestScore:
                yvals.append(localBestScore)
                xvals.append(time.time()-startTime)
        
        print "trajectory number:",i+1,"out of",numTraj
        done = False
        while not done:
            allCurrGates = graph.getAllGates()
            #Only consider swaps that would actually change the graph at 
            #least one of the genes in a pair must be in the graph already
            swappablePairs = []
            for pair in standardAllRepressorSwappingPairs:
                if pair[0] in allCurrGates or pair[1] in allCurrGates:
                    swappablePairs.append(pair)
            random.shuffle(swappablePairs)
            for pair in swappablePairs:
                #print "On pair",swappablePairs.index(pair)+1,"out of",len(swappablePairs)
                graph.swapGates(pair[0],pair[1])
                
                #make sure the we didnt introduce a repeated gene into the graph
                if graph.hasRepeatedRepressor():
                    graph.swapGates(pair[0],pair[1])
                    continue
                graph.writeToJson(scratchFile)
                tempScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
                tempResults = getScoreFromDict(outputNames, tempScoreDict, smallestScoreAllowed)
                #if there was an improvement that passed the criteria for being
                #a valid assignment save the new score and graph.
                if tempResults[1] and tempResults[0]>localBestScore:
                    localBestGraph = copy.deepcopy(graph)
                    localBestScore = tempResults[0]
                    localBestScoreDict = tempScoreDict
                    if localBestScore>globalBestScore:
                        yvals.append(localBestScore)
                        xvals.append(time.time()-startTime)
                    
                #if there was an improvement then break the loop whether or
                #not it passes the criteria for being a valid assignment
                if tempResults[0]>currScore:
                    print "improved"
                    currScore = tempResults[0]
                    break
                #if there was not an improvement switch back
                elif tempResults[0]<=currScore:
                    graph.swapGates(pair[0],pair[1])
                #we are done with trajectory if we make it to the end of the list
                if swappablePairs.index(pair)==len(swappablePairs)-1:
                    done = True
                    
        #check if the local trajectory was better than the global one
        if localBestGraph!=None and localBestScore>globalBestScore:
            globalBestGraph = localBestGraph
            globalBestScore = localBestScore
            globalBestScoreDict = localBestScoreDict
            
#    print globalBestScore
#    print globalBestGraph
    xvals.append((time.time()-startTime))
    yvals.append(globalBestScore)
    endTime = time.time()
    print "This took",(endTime-startTime),"seconds"
#    print globalBestGraph
    graphResults(xvals,yvals)    
    
    return globalBestGraph,globalBestScore,globalBestScoreDict

def findOptimalAssortmentHillTimed(netlist,Libraries=Libraries, smallestScoreAllowed=10,numTraj=5,maxTime=300):
    xvals = [0]
    yvals = [0]
    startTime = time.time()
    #random.seed(0)
    scratchFile = "JsonFiles/optimizingScratchWorkHillTimed.json"
    graph, allGates, Inputs, Repressors, Outputs = OptimalCircuit.makeGraphFromNetlist(netlist, True)
    if graph.hasLoop():
        print "This circuit is sequential and cannot be scored."
        return
    allInputs = makeInputs(Libraries[0])
    allRepressors = makeRepressors(Libraries[1])
    allOutputs = makeOutputs(Libraries[2])
    numInputs = len(Inputs)
    numRepressors = len(Repressors)
    numOutputs = len(Outputs)
    chosenInputs = selectInputs(allInputs,numInputs)
    chosenOutputs = selectOutputs(allOutputs,numOutputs)
    for i in range(numInputs):
        graph.swapGates(chosenInputs[i],Inputs[i])
    outputNames = []
    for i in range(numOutputs):
        graph.swapGates(chosenOutputs[i], Outputs[i])
        outputNames.append(chosenOutputs[i].getName())

    #initialize the global best
    globalBestGraph = None
    globalBestScore = 0
    globalBestScoreDict = None
    #get all possible swaps
    standardAllRepressorSwappingPairs = list(itertools.combinations(allRepressors,2))
    
    #start the trajectory
    for i in range(numTraj):
        if (time.time()-startTime >= maxTime):
            break
        #randomly assign the starting point
        chosenStartRepressors = selectStartingRepressors(allRepressors,numRepressors)
        for j in range(numRepressors):
            graph.swapGates(chosenStartRepressors[j], Repressors[j])
        Repressors = chosenStartRepressors
        while graph.hasRepeatedRepressor():
            chosenStartRepressors = selectStartingRepressors(allRepressors,numRepressors)
            for j in range(numRepressors):
                graph.swapGates(chosenStartRepressors[j], Repressors[j])
            Repressors = chosenStartRepressors

        #initialize the local best
        localBestGraph = None
        localBestScore = 0
        localBestScoreDict = None
        #get the score of this starting point
        graph.writeToJson(scratchFile)
        firstScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
        firstResults = getScoreFromDict(outputNames, firstScoreDict, smallestScoreAllowed)
        currScore = firstResults[0]
        if firstResults[1]:
            localBestGraph = copy.deepcopy(graph)
            localBestScore = firstResults[0]
            localBestScoreDict = firstScoreDict
            if localBestScore>globalBestScore:
                yvals.append(localBestScore)
                xvals.append(time.time()-startTime)
        
        print "trajectory number:",i+1,"out of",numTraj
        done = False
        while not done:
            allCurrGates = graph.getAllGates()
            #Only consider swaps that would actually change the graph
            swappablePairs = []
            for pair in standardAllRepressorSwappingPairs:
                if pair[0] in allCurrGates or pair[1] in allCurrGates:
                    swappablePairs.append(pair)
            random.shuffle(swappablePairs)
            for pair in swappablePairs:
                if (time.time()-startTime >= maxTime):
                    done = True
                    break
                #print "On pair",swappablePairs.index(pair)+1,"out of",len(swappablePairs)
                graph.swapGates(pair[0],pair[1])
                if graph.hasRepeatedRepressor():
                    graph.swapGates(pair[0],pair[1])
                    continue
                graph.writeToJson(scratchFile)
                tempScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
                tempResults = getScoreFromDict(outputNames, tempScoreDict, smallestScoreAllowed)
                #if there was an improvement that passed the criteria for being
                #a valid assignment save the new score and graph.
                if tempResults[1] and tempResults[0]>localBestScore:
                    localBestGraph = copy.deepcopy(graph)
                    localBestScore = tempResults[0]
                    localBestScoreDict = tempScoreDict
                    if localBestScore>globalBestScore:
                        yvals.append(localBestScore)
                        xvals.append(time.time()-startTime)
                    
                #if there was an improvement then break the loop whether or
                #not it passes the criteria for being a valid assignment
                if tempResults[0]>currScore:
                    print "improved"
                    currScore = tempResults[0]
                    break
                #if there was not an improvement switch back
                elif tempResults[0]<=currScore:
                    graph.swapGates(pair[0],pair[1])
                #we are done with trajectory if we make it to the end of the list
                if swappablePairs.index(pair)==len(swappablePairs)-1:
                    done = True
        #check if the local trajectory was better than the global one
        if localBestGraph!=None and localBestScore>globalBestScore:
            globalBestGraph = localBestGraph
            globalBestScore = localBestScore
            globalBestScoreDict = localBestScoreDict

#    print globalBestScore
#    print globalBestGraph
    xvals.append((time.time()-startTime))
    yvals.append(globalBestScore)
    endTime = time.time()
    print "This took",(endTime-startTime),"seconds"
    graphResults(xvals,yvals)
    return globalBestGraph,globalBestScore,globalBestScoreDict,xvals,yvals

def findOptimalAssortmentRandTimed(netlist,Libraries=Libraries, smallestScoreAllowed=10,maxTime=300):
    startTime = time.time()
    #random.seed(0)
    scratchFile = "JsonFiles/optimizingScratchWorkRandTimed.json"
    graph, allGates, Inputs, Repressors, Outputs = OptimalCircuit.makeGraphFromNetlist(netlist, True)
    if graph.hasLoop():
        print "This circuit is sequential and cannot be scored."
        return
    allInputs = makeInputs(Libraries[0])
    allRepressors = makeRepressors(Libraries[1])
    allOutputs = makeOutputs(Libraries[2])
    numInputs = len(Inputs)
    numRepressors = len(Repressors)
    numOutputs = len(Outputs)
    chosenInputs = selectInputs(allInputs,numInputs)
    chosenOutputs = selectOutputs(allOutputs,numOutputs)
    for i in range(numInputs):
        graph.swapGates(chosenInputs[i],Inputs[i])
    outputNames = []
    for i in range(numOutputs):
        graph.swapGates(chosenOutputs[i], Outputs[i])
        outputNames.append(chosenOutputs[i].getName())

    #initialize the global best
    globalBestGraph = None
    globalBestScore = 0
    globalBestScoreDict = None
    xvals = [0]
    yvals = [0]
    #get all possible swaps
    while (time.time()-startTime)<maxTime:
        chosenRepressors = selectStartingRepressors(allRepressors,numRepressors)
        for j in range(numRepressors):
            graph.swapGates(chosenRepressors[j], Repressors[j])
        Repressors = chosenRepressors
        
        #fix cases where we have same repressor
        while graph.hasRepeatedRepressor():
            chosenRepressors = selectStartingRepressors(allRepressors,numRepressors)
            for j in range(numRepressors):
                graph.swapGates(chosenRepressors[j], Repressors[j])
            Repressors = chosenRepressors
        
        graph.writeToJson(scratchFile)
        tempScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
        tempResults = getScoreFromDict(outputNames, tempScoreDict, smallestScoreAllowed)
        
        #check if this circuit is better than the best found so far
        if tempResults[1] and tempResults[0]>=globalBestScore:
            globalBestGraph = copy.deepcopy(graph)
            globalBestScore = tempResults[0]
            globalBestScoreDict = tempScoreDict
            xvals.append((time.time()-startTime))
            yvals.append(globalBestScore)
    xvals.append((time.time()-startTime))
    yvals.append(globalBestScore)
    
    print globalBestGraph
    graphResults(xvals,yvals)
    return globalBestGraph,globalBestScore,globalBestScoreDict,xvals,yvals
    
def compareHillClimbToRandom(netlist,Libraries=Libraries, smallestScoreAllowed=3,numTraj=5,maxTime=300):
    """
    Compares the score vs time progression for the random and the hill climbing algorithm.
    """
    a = findOptimalAssortmentHillTimed(netlist,Libraries, smallestScoreAllowed=smallestScoreAllowed,numTraj=numTraj,maxTime=maxTime)
    b = findOptimalAssortmentRandTimed(netlist,Libraries, smallestScoreAllowed=smallestScoreAllowed,maxTime=maxTime)
    x1vals = a[3]
    y1vals = a[4]
    x2vals = b[3]
    y2vals = b[4]
    plt.figure()
    ymax = max([max(y1vals),max(y2vals)])*1.1
    xmax = max([max(x1vals),max(x2vals)])*1.1
    plt.plot(x1vals,y1vals,color="blue",linestyle="-",marker="o",label="HillClimbing")
    plt.plot(x2vals,y2vals,color="red",linestyle="-",marker="o",label="Random")
    plt.xlabel('Time (sec)')
    plt.ylabel('Score')
    plt.title('Score Progression')
    plt.axis(xmin=0, xmax=xmax, ymin=0, ymax = ymax)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
#--------Helper Functions --------
#Make gates from libraries takes in file directories to the libraries and 
#generates a list of gates from the information in the files.
def makeInputs(inputLibrary):
    """
    returns a list of inputs gates generated from the information in the 
    inputs library json file
    """
    myFile = open(inputLibrary,'r')
    inputsFromFile = json.load(myFile)
    myFile.close()
    
    allInputs = []
    for gateInfo in inputsFromFile:
#        Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
        tempGate = Gate.Gate(str(gateInfo["NAME"]),"Input",gateInfo["Km"],gateInfo["n"],gateInfo["Pmin"],gateInfo["Pmax"],None,None)
        allInputs.append(tempGate)
    return allInputs
    
def makeOutputs(outputLibrary):
    """
    returns a list of outputs gates generated from the information in the 
    outputs library json file
    """
    myFile = open(outputLibrary,'r')
    outputsFromFile = json.load(myFile)
    myFile.close()
    
    allOutputs = []
    for gateInfo in outputsFromFile:
#        Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
        tempGate = Gate.Gate(str(gateInfo["NAME"]),"Output",None,None,None,None,gateInfo["halfLife"],None)
        allOutputs.append(tempGate)
    return allOutputs

def makeRepressors(repressorLibrary):
    """
    returns a list of repressor gates generated from the information in the 
    repressor library json file
    """
    myFile = open(repressorLibrary,'r')
    repressorsFromFile = json.load(myFile)
    myFile.close()
    
    allRepressors = []
    
    for gateInfo in repressorsFromFile:
#        Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
        tempGate = Gate.Gate(str(gateInfo["NAME"]),"Repressor",gateInfo["Km"],gateInfo["n"],gateInfo["Pmin"],gateInfo["Pmax"],gateInfo["halfLife"],None)
        allRepressors.append(tempGate)
    return allRepressors
    
#Choose the genes that will be used in the circuit.
def selectInputs(inputsList, numInputs):
    """
    Takes in some integer and returns a list with that many input gates.
    """    
    #TODO
    #return inputsList[:numInputs]
    return random.sample(inputsList,numInputs)

def selectOutputs(outputsList,numOutputs):
    """
    Takes in some integer and returns a list with that many output gates.
    """    
    #TODO
    #return outputsList[:numOutputs]
    return random.sample(outputsList,numOutputs)
    
def selectStartingRepressors(repressorsList,numRepressors):
    """
    Takes in some integer and returns a list with that many repressors gates.
    """    
    #TODO
    return random.sample(repressorsList,numRepressors)
    
def graphResults(xvals,yvals):
    """
    given a set of points where the xvals are time points and the yvals are 
    scores, graphs the score progression.
    """
    plt.figure()
    ymax = max(yvals)*1.1
    xmax = max(xvals)*1.1
    plt.plot(xvals,yvals,color="blue",linestyle="-",marker="o")
    plt.xlabel('Time (sec)')
    plt.ylabel('Score')
    plt.title('Score Progression')
    plt.axis(xmin=0, xmax=xmax, ymin=0, ymax = ymax)
    plt.show()
    

def getScoreFromDict(outputNames, scoreDict, smallestScoreAllowed=10):
    """
    given a dictionary of scores for each gate and the names of the outputs and
    the minimum allowed score, returns the score of the output followed by true
    if all intermediate gates are greater than the minimum allowed score and
    false if not.
    """
    score = 0
    isAcceptable = True
    for key in scoreDict:
        if scoreDict[key]<smallestScoreAllowed:
            isAcceptable = False
        if key in outputNames:
            score += scoreDict[key]
    
    return (score,isAcceptable)