# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:17:06 2015

@author: Arinze
"""
import json
import Gate
import Graph
import itertools
import OptimalCircuit
import General
import copy
import random
inputsDir = "JsonFiles/Libraries/InputLibrary2.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary2.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary2.json"
Libraries = [inputsDir, repressorsDir, outputsDir]
truthValueExampleFileLoc = "JsonFiles/SplitByTruthValue_OR/01101001.json" #01101001
placeToSaveTruthValueExample = "JsonFiles/01101001Graph.json"
exampleNetlist = ['NOR(W1,a,b)','NOR(W2,W1,c)','NOR(W3,a,c,d)','NOR(W4,W1,a,d)','NOR(W5,W2,W3,W4)']

def findOptimalAssortmentHill(netlist,Libraries, smallestScoreAllowed=100,numTraj=20):
    random.seed(0)
    scratchFile = "JsonFiles/optimizingScratchWork.json"
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
    #get all possible swaps
    standardAllRepressorSwappingPairs = list(itertools.combinations(allRepressors,2))
    
    #start the trajectory
    for i in range(numTraj):
        #randomly assign the starting point
        chosenStartRepressors = selectStartingRepressors(allRepressors,numRepressors)
        for j in range(numRepressors):
            graph.swapGates(chosenStartRepressors[j], Repressors[j])
        Repressors = chosenStartRepressors

        #initialize the local best
        localBestGraph = None
        localBestScore = 0
        #get the score of this starting point
        graph.writeToJson(scratchFile)
        firstScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
        firstResults = getScoreFromDict(outputNames, firstScoreDict, smallestScoreAllowed)
        currScore = firstResults[0]
        if firstResults[1]:
            localBestGraph = copy.deepcopy(graph)
            localBestScore = firstResults[0]
        
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
                graph.swapGates(pair[0],pair[1])
                graph.writeToJson(scratchFile)
                tempScoreDict = General.generateDynamicCircuitGraphs(scratchFile, False, False, False)
                tempResults = getScoreFromDict(outputNames, tempScoreDict, smallestScoreAllowed)
                #if there was an improvement that passed the criteria for being
                #a valid assignment save the new score and graph.
                if tempResults[1] and tempResults[0]>localBestScore:
                    localBestGraph = copy.deepcopy(graph)
                    localBestScore = tempResults[0]
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
#            print "this traj was better"
            globalBestGraph = localBestGraph
            globalBestScore = localBestScore
#        else:
#            print "this traj was not better"
    print globalBestScore
    print globalBestGraph
    return globalBestGraph,globalBestScore



def getScoreFromDict(outputNames, scoreDict, smallestScoreAllowed=100):
    score = 0
    isAcceptable = True
    for key in scoreDict:
        if scoreDict[key]<smallestScoreAllowed:
            isAcceptable = False
        if key in outputNames:
            score += scoreDict[key]
    
    return (score,isAcceptable)
#--------Helper Functions --------
def selectInputs(inputsList, numInputs):
    """
    Takes in some integer and returns a list with that many input gates.
    """    
    #TODO
    return inputsList[:numInputs]

def selectOutputs(outputsList,numOutputs):
    """
    Takes in some integer and returns a list with that many output gates.
    """    
    #TODO
    return outputsList[:numOutputs]
def selectStartingRepressors(repressorsList,numRepressors):
    """
    Takes in some integer and returns a list with that many repressors gates.
    """    
    #TODO
    return random.sample(repressorsList,numRepressors)
#Make gates from libraries takes in file directories to the libraries and 
#generates a list of gates from the information in the files.
def makeInputs(inputLibrary):
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
    myFile = open(outputLibrary,'r')
    outputsFromFile = json.load(myFile)
    myFile.close()
    
    allOutputs = []
    for gateInfo in outputsFromFile:
#        Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
        tempGate = Gate.Gate(str(gateInfo["NAME"]),"Repressor",gateInfo["Km"],gateInfo["n"],gateInfo["Pmin"],gateInfo["Pmax"],gateInfo["halfLife"],None)
        allOutputs.append(tempGate)
    return allOutputs
    
def makeRepressors(repressorLibrary):
    myFile = open(repressorLibrary,'r')
    repressorsFromFile = json.load(myFile)
    myFile.close()
    
    allRepressors = []
    
    for gateInfo in repressorsFromFile:
#        Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
        tempGate = Gate.Gate(str(gateInfo["NAME"]),"Output",None,None,None,None,gateInfo["halfLife"],None)
        allRepressors.append(tempGate)
    return allRepressors