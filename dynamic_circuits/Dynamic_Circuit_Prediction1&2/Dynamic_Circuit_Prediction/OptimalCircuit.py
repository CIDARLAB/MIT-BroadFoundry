# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 12:09:16 2015

@author: Arinze
"""

import Graph
import Gate
import Wire
import inputs
import json
import General
import GeneralJsonIO
import inputs
import update
import re
import itertools
import time

inputsDir = "JsonFiles/Libraries/InputLibrary2.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary2.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary2.json"
Libraries = [inputsDir, repressorsDir, outputsDir]
    
repressilatorFileLoc = "JsonFiles/repressilator.json"
placeToSaveRepressilator = "JsonFiles/repressilatorGraph.json"
#
truthValueExampleFileLoc = "JsonFiles/SplitByTruthValue_OR/01101001.json" #01101001
placeToSaveTruthValueExample = "JsonFiles/01101001Graph.json"
netlistLoc = "JsonFiles/standardOR.json"
intermediateFile = "JsonFiles/intermediateGraph.json"
#Latches
Gated_D_Latch_FileLoc = "JsonFiles/SR_Latches/Gated_D_Latch.json"
Gated_D_Latch_placeToSave = "JsonFiles/SR_Latches/Gated_D_Latch_Graph.json"
Gated_D_Latch_noCInv_FileLoc = "JsonFiles/SR_Latches/Gated_D_Latch_noCInv.json"
Gated_D_Latch_noCInv_placeToSave = "JsonFiles/SR_Latches/Gated_D_Latch_noCInv_Graph.json"
Negative_Edge_Triggered_D_Flip_Flop_FileLoc = "JsonFiles/SR_Latches/Negative_Edge_Triggered_D_Flip_Flop.json"
Negative_Edge_Triggered_D_Flip_Flop_placeToSave = "JsonFiles/SR_Latches/Negative_Edge_Triggered_D_Flip_Flop_Graph.json"
Transparent_D_Latch_noDinv_FileLoc = "JsonFiles/SR_Latches/Transparent_D_Latch_noDinv.json"
Transparent_D_Latch_noDinv_placeToSave = "JsonFiles/SR_Latches/Transparent_D_Latch_noDinv_Graph.json"
Transparent_D_Latch_noInpinv_FileLoc = "JsonFiles/SR_Latches/Transparent_D_Latch_noInpinv.json"
Transparent_D_Latch_noInpinv_placeToSave = "JsonFiles/SR_Latches/Transparent_D_Latch_noInpinv_Graph.json"
#TestCase
Test1_FileLoc = "JsonFiles/SR_Latches/Test1.json"
Test1_placeToSave = "JsonFiles/SR_Latches/Test1_Graph.json"
Test2_FileLoc = "JsonFiles/SR_Latches/Test2.json"
Test2_placeToSave = "JsonFiles/SR_Latches/Test2_Graph.json"
#
truthValue = "00010001"
inputNames = ["IPTG1","IPTG2","IPTG3"]
outputNames = ["FP1"]

def makeGraphFromNetlist(inputNetlist, useDefaultInput):
    """
    Takes in a circuit in the string form like '((a.b).c)' and returns a graph 
    created from that string, a list of all the gates in that graph, lists of the
    inputs, repressors, and outputs individually from all the gates.
    This assumes all inputs are activators and allows for intermediate genes to
    be activators. Treats all outputs are buffers or ORs.
    """
    if type(inputNetlist) == str:
        myFile = open(inputNetlist,'r')
        netlist = json.load(myFile)
        myFile.close()
    elif type(inputNetlist) == list:
        netlist = inputNetlist
    
    #Initialize the lists of components
    Inputs = []
    Intermediates = []
    Outputs = []
    allWires = []
    fanOutWireToGateDict = {}
    allGatesParsed = {}
    allFanOutNames = []
    allFanInNames = []
    allWireNames = []
    
    for gateInfo in netlist:
        gateInfo = gateInfo.replace(")","")
        gateInfo = gateInfo.replace("(",",")
        temp = gateInfo.split(",")
        allGatesParsed[temp[1]] = []
        allFanOutNames.append(temp[1])
        if temp[1] not in allWireNames:
            allWireNames.append(temp[1])
        for i in range(len(temp)):
            if i!=1:
                allGatesParsed[temp[1]].append(temp[i])
            if i>1:
                if temp[i] not in allFanInNames:
                    allFanInNames.append(temp[i])
                if temp[i] not in allWireNames:
                    allWireNames.append(temp[i])
             
    givenInputs = []
    givenOutputs = []
    givenIntermediates = []
    for i in allFanInNames:
        if i not in allFanOutNames:
            givenInputs.append(i)
    givenInputs.sort()
    #givenInputs.reverse()
    for i in allFanOutNames:
        if i not in allFanInNames:
            givenOutputs.append(i)
    for i in allWireNames:
        if (i not in givenInputs) and (i not in givenOutputs):
            givenIntermediates.append(i)
    #Make the inputs
    numInputs = len(givenInputs)
    initPeriod = 500.0*2**(numInputs)
    period = initPeriod
    inputAliases = {}
    if useDefaultInput:
        inputBinaries = makeInputBin(numInputs)
    for i in range(numInputs):
        standardGateName = "IN"+str(i+1)
        #Create the input gate and set the equation to use and its expected protein
        #and promoter truth values which should be the same for an input.
        #Gate(name,Km,n,gateType,mB,pB,m_halfLife, p_halfLife)
        tempInput = Gate.Gate(standardGateName,1,2,'Input',30,None,None,None)
        if useDefaultInput:
            tempTV = ''.join(inputBinaries[i])
            tempInput.setInputType(['inputs.specInput',100.0,0.0,initPeriod,tempTV])
        else:
            tempTV = raw_input(givenInputs[i]+"'s input truthvalue: ")
            tempInput.setInputType(['inputs.squInput',period,10,period/2.0,0])
        
        
        tempInput.setExpectedProtein(list(tempTV))
        inputAliases[givenInputs[i]]=standardGateName
        #So the period of the next input alternates twice as frequently as this
        period /= 2.0
        #Add it to the list of inputs. Note the name of the fanOut wire as being
        #from that gate and replace all occurences of that letter with the wire
        #name associated with that gate.
        Inputs.append(tempInput)
        fanOutWireToGateDict[givenInputs[i]] = tempInput
    
    #Make the intermediate gates
    gateCount = 1
    for i in range(len(givenIntermediates)):
        standardGateName = "G" + str(gateCount)
        givenWireName = givenIntermediates[i]
        if allGatesParsed[givenWireName][0]=="NOR" or allGatesParsed[givenWireName][0]=="NOT": 
            tempIntermediate = Gate.Gate(standardGateName,1000,2,'Repressor',30,20,2,10)
        elif allGatesParsed[givenWireName][0]=="OR" or allGatesParsed[givenWireName][0]=="BUF": 
            tempIntermediate = Gate.Gate(standardGateName,1000,2,'Activator',30,20,2,10)
        Intermediates.append(tempIntermediate)
        fanOutWireToGateDict[givenWireName] = tempIntermediate
        gateCount += 1
    
    #Make outputs
    numOutputs = 1
    outputAliases = {}
    for i in range(len(givenOutputs)):
        givenWireName = givenOutputs[i]
        standardOutName = "Y"+str(numOutputs)
        if givenWireName==standardOutName:
            standardOutName = "Y"+str(numOutputs)+"r"
        if allGatesParsed[givenWireName][0]=="NOR" or allGatesParsed[givenWireName][0]=="NOT": 
            standardGateName = "G" + str(gateCount)
            tempIntermediate = Gate.Gate(standardGateName,1000,2,'Repressor',30,20,2,10)
            Intermediates.append(tempIntermediate)
            fanOutWireToGateDict[givenWireName] = tempIntermediate
            allGatesParsed[standardOutName] = ["BUF",givenWireName]      
            
            tempOutput = Gate.Gate(standardOutName,None,None,'Output',None,20,2,10)
            Outputs.append(tempOutput)
            fanOutWireToGateDict[standardOutName] = tempOutput
            outputAliases[givenWireName] = standardOutName
            gateCount += 1

        elif allGatesParsed[givenWireName][0]=="OR" or allGatesParsed[givenWireName][0]=="BUF": 
            tempOutput = Gate.Gate(standardOutName,None,None,'Output',None,20,2,10)
            Outputs.append(tempOutput)
            fanOutWireToGateDict[givenWireName] = tempOutput
            outputAliases[givenWireName] = standardOutName
        else:
            print "Something isn't right."
        numOutputs += 1

    #Make and connect wires
    wireCount = 1
    wireAliases = {}
    for givenWireName in allGatesParsed.keys():
        #get the gate associated with that wire
        tempGate = fanOutWireToGateDict[givenWireName]
        #Ignore the first element which tells the type of gate
        fanInWireNames = allGatesParsed[givenWireName][1:]
        for wireName in fanInWireNames:
            if wireName not in wireAliases:
                wireAliases[wireName] = "W"+str(wireCount)
                wireCount += 1
            wire = Wire.Wire(wireAliases[wireName])
            allWires.append(wire)
            associatedFromGate = fanOutWireToGateDict[wireName]
            tempGate.addFanInWire(wire)
            associatedFromGate.addFanOutWire(wire)
    allGates = Inputs + Intermediates + Outputs
    
    if numInputs > 0:
        graph = Graph.Graph(0,initPeriod,10*initPeriod)
    else:
        #arbitrary
        graph = Graph.Graph(0,1000,10000)
    for wire in allWires:
        graph.addWire(wire)
    for gate in allGates:
        graph.addGate(gate)
    print inputAliases
    print outputAliases
    return graph, allGates, Inputs, Intermediates, Outputs

    
def performSwaps(graph,allGates,Inputs,allGatesr,Inputsr):
    """
    Given a graph, a set of gates in the graph and a new set of gates, it will swap
    in the new gates for the old gates. The number of each type of gate must be
    the same.
    """
    
    #Make swaps in graph
#    print "Performing swaps"
    for i in range(len(allGates)):
        if allGatesr[i].getGateType() != "Output":
            allGatesr[i].setFanOut(allGates[i].getFanOut())
        if allGatesr[i].getGateType() != "Input":
            allGatesr[i].setFanIn(allGates[i].getFanIn())
    for i in range(len(Inputsr)):
        Inputsr[i].setInputType(Inputs[i].getInputType())
        Inputsr[i].setExpectedProtein(Inputs[i].getExpectedProtein())
    #add new gates to graph 
    for gate in allGatesr:
        graph.addGate(gate)
    #remove old gates from graph
    for gate in allGates:
        graph.removeGate(gate)

def wrapperForNetlist(fileLoc, placeToSave, makeBarGraph=True,makeOtherGraphs=True,useDefaultInput=True):
    """
    Takes in a netlist of a circuit, and a directory to store the intermediate json file.
    Prints the graphs generated for mRNA and protein concentrations generated
    by General. Returns the score for the circuit with the assigned gates.
    """
    #Make Graph from string
    graph, allGates, Inputs, Intermediates, Outputs = makeGraphFromNetlist(fileLoc,useDefaultInput)

    isSequential = False
    print graph
    if graph.hasLoop():
        isSequential = True
        #So things do not start at equilibrium.
        Intermediates[0].setInitialProtein(1.0)
    #Write to Json File
    graph.writeToJson(placeToSave)
    #Make graphs from JsonFile
    scores = General.generateDynamicCircuitGraphs(placeToSave, makeBarGraph, makeOtherGraphs, isSequential)
    return scores
    
         
def makeInputBin(numInputs):
    binSize = 2**numInputs
    base = "01"
    allInputBin = []
    for i in range(numInputs):
        
        tempInput = base
        while len(tempInput)<binSize:
            tempInput += base
        nextLen = len(base)*2
        while len(base)<nextLen:
            base = "0" + base + "1"
        allInputBin.append(list(tempInput))
        #allInputBin.append(tempInput)
    allInputBin.reverse()
    return allInputBin


def findBestForTruthValue(truthValue, inputNames, outputNames, intermediateFile = intermediateFile, netlistLoc = netlistLoc, Libraries = Libraries, minAllowedScore = 100, makeBarGraph = True):
    """
    Assumes the length of the inputNames is the same as the number of inputs 
    necessary and the same for outputs. Returns the best circuit repressors and
    netlist for the desired truth value.
    """
    startTime = time.time()
    myFile = open(netlistLoc,'r')
    data = json.load(myFile)
    myFile.close()
    
    bestNetlist = None
    bestRepressorOrder = None
    bestIntermediateGateNames = None
    bestScore = 0
    bestScoreForNetlist = bestScore
    
    allNetlists = data[truthValue]["netlists"]
    initInputsr = makeInputs(Libraries[0],inputNames)
    initOutputsr = makeOutput(Libraries[2],outputNames)
    for netlist in allNetlists:
        print "best so far:", bestScoreForNetlist
        if bestScoreForNetlist!=0 and bestScoreForNetlist<bestScore:
            print "Beat last netlist"
        bestScoreForNetlist = bestScore
        graph, allGates, Inputs, Intermediates, Outputs = makeGraphFromNetlist(netlist,True)
        print graph
        IntermediatesrList = makeRepressorsGroups(Libraries[1],len(Intermediates))
        for IntermediatesrTuple in IntermediatesrList:
            #-------------------------------
            graph, allGates, Inputs, Intermediates, Outputs = makeGraphFromNetlist(netlist,True)
            #This will tell how the repressor name order corresponds to gate number
            #from the netlist
            IntermediateGateNames = []
            for gate in Intermediates:
                IntermediateGateNames.append(gate.getName())
            #These two lines should not affect anything if the correct number of
            #names are specified        
            Inputsr = initInputsr[0:len(Inputs)]
            Outputsr = initOutputsr[0:len(Outputs)]
            #---------------------------------
            Intermediatesr = list(IntermediatesrTuple)
            IntermediatesrNames = []
            for gate in Intermediatesr:
                IntermediatesrNames.append(gate.getName())
            print IntermediatesrNames
            allGatesr = Inputsr + Intermediatesr + Outputsr
            performSwaps(graph, allGates, Inputs, allGatesr, Inputsr)
#            print graph
            graph.writeToJson(intermediateFile)
            if graph.hasLoop():
                Intermediatesr[0].setInitialProtein(1.0)
                makeBarGraph = False
                print "Problem"
                
            score = General.generateDynamicCircuitGraphs(intermediateFile, makeBarGraph)
            if score == None:
                print "Problem"
                score = {None:-1}
            allGates = score.keys()
            isGood = True
            for gate in allGates:
                if score[gate]<minAllowedScore:
                    isGood = False
                    break
            totScore = 0
            if isGood:
                for name in outputNames:
                    totScore += score[name]
            if totScore>bestScore:
                bestScore = totScore
                bestNetlist = netlist
                bestRepressorOrder = IntermediatesrNames
                bestIntermediateGateNames = IntermediateGateNames

#            print totScore
#            print bestScore
#            pause = raw_input("pause")
    gateAliases = {}
    for i in range(len(bestIntermediateGateNames)):
        gateAliases[bestIntermediateGateNames[i]] = bestRepressorOrder[i]
    endTime = time.time()
    print "This took",round(endTime-startTime,2),"seconds."
    return bestNetlist, gateAliases, bestScore
            

#Make gates from libraries
def makeInputs(inputLibrary,inputNames):
    myFile = open(inputLibrary,'r')
    inputsFromFile = json.load(myFile)
    myFile.close()
    
    inputsFromFileNames = []
    for gateInfo in inputsFromFile:
        inputsFromFileNames.append(gateInfo["NAME"])
    
    Inputsr = []
    for name in inputNames:
        inputIndex = inputsFromFileNames.index(name)
        gateInfo = inputsFromFile[inputIndex]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),gateInfo["Km"],gateInfo["n"],"Input",gateInfo["mB"],None,None,None)
        Inputsr.append(tempGate)
    return Inputsr
    
def makeOutput(outputLibrary,outputNames):
    myFile = open(outputLibrary,'r')
    outputsFromFile = json.load(myFile)
    myFile.close()
    
    outputsFromFileNames = []
    for gateInfo in outputsFromFile:
        outputsFromFileNames.append(gateInfo["NAME"])
        
    Outputsr = []
    for name in outputNames:
        outputIndex = outputsFromFileNames.index(name)
        gateInfo = outputsFromFile[outputIndex]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),None,None,"Output",None,gateInfo["pB"],gateInfo["m_halfLife"],gateInfo["p_halfLife"])
        Outputsr.append(tempGate)
    return Outputsr
    
def makeRepressorsGroups(repressorLibrary,numRepressors):
    myFile = open(repressorLibrary,'r')
    repressorsFromFile = json.load(myFile)
    myFile.close()
    
    allRepressors = []
    
    for gateInfo in repressorsFromFile:
        tempGate = Gate.Gate(str(gateInfo["NAME"]),gateInfo["Km"],gateInfo["n"],"Repressor",gateInfo["mB"],gateInfo["pB"],gateInfo["m_halfLife"],gateInfo["p_halfLife"])
        allRepressors.append(tempGate)
    #Make list of permutations
    IntermediatesrList = list(itertools.permutations(allRepressors, numRepressors))
    return IntermediatesrList

def findSmallestParentheses(circuit):
    """
        Returns the indeces of a pair of parentheses that do not contain
        any parentheses in it
    """
    L = len(circuit)
    for i in xrange(L):
        if circuit[i]=="(":
            startIndex = i
        if circuit[i]==")":
            endIndex = i
            try:
                return (startIndex,endIndex)
            except NameError:
                print "parentheses are not properly matched"
                return None
    return None
    
def convertCircuitStringToNetlist(circuitString):
    """
    Converts the string representation of a circuit to a netlist.
    """
    operators = ["~", "&", "@","+","^",".","="]
    operatorNames = {"~":"NOT", "&":"AND", "@":"NAND","+":"OR","^":"XOR",".":"NOR","=":"XNOR"}
    
    netlist = []
    wireCount = 1
    indeces = findSmallestParentheses(circuitString)
    while indeces!=None:
        wireName = "W"+str(wireCount)
        subCircuit = circuitString[indeces[0]:indeces[1]+1]
        circuitString = circuitString.replace(subCircuit,wireName)
        
        subCircuit = subCircuit[1:-1]
        opUsed = "x"
        for op in operators:
            components = subCircuit.split(op)
            if len(components)>1:
                opUsed = op
                break
        if opUsed == "x":
            "Error"
            assert False
        if(components[0] == ""):
            components.pop(0)
        term = operatorNames[opUsed]+"("+wireName+","
        for c in components:
            term = term + c + ","
        term = term[:-1] + ")"
        
        netlist.append(term)
        indeces = findSmallestParentheses(circuitString)
        wireCount += 1
    return netlist
    