# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 12:09:16 2015

@author: Arinze
"""

import Graph
import Gate
import Wire
import json
import General
import Optimize

inputsDir = "JsonFiles/Libraries/InputLibrary.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary.json"
Libraries = [inputsDir, repressorsDir, outputsDir]
    
repressilatorFileLoc = "JsonFiles/repressilator.json"
placeToSaveRepressilator = "JsonFiles/repressilatorGraph.json"
repressilatorFileLoc2 = "JsonFiles/repressilator2.json"
placeToSaveRepressilator2 = "JsonFiles/repressilatorGraph2.json"

netlistLoc = "JsonFiles/standardOR.json"
intermediateFile = "JsonFiles/intermediateGraph.json"

truthValue = "00010001"
inputNames = ["IPTG1","IPTG2","IPTG3"]
outputNames = ["FP1"]

#useDefaultInput creates a combinational truth table.
#might not want to use default for sequential, a user-defined waveform might be desired instead.
def makeGraphFromNetlist(inputNetlist, useDefaultInput=True):
    """
    Takes in a circuit in the string form like '((a.b).c)' and returns a graph 
    created from that string, a list of all the gates in that graph, lists of the
    inputs, repressors, and outputs individually from all the gates.
    This assumes all inputs are activators and allows for intermediate genes to
    be activators. Treats all outputs are buffers or ORs.
    """
    if type(inputNetlist) == str: #file path, read netlist from file
        myFile = open(inputNetlist,'r')
        netlist = json.load(myFile)
        myFile.close()
    elif type(inputNetlist) == list: #netlis = list of strings
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
    
    #Parses through netlist to obtain gate attributes; netlist is a list 
    #composed of information about the gateType(fanOut, fanIn [there can be 
    #more than 1]}
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
    #Inputs, gates with fanIn wires with no FanOut origin
    for i in allFanInNames:
        if i not in allFanOutNames:
            givenInputs.append(i)
    givenInputs.sort()
    #Outputs, gates with fanOut wires with no FanIn destination
    for i in allFanOutNames:
        if i not in allFanInNames:
            givenOutputs.append(i)
    #Intermediate nodes, middle gates with fanIn and fanOut wires
    for i in allWireNames:
        if (i not in givenInputs) and (i not in givenOutputs):
            givenIntermediates.append(i)
            
    #Creates the logic-input function as a set of square waves
    numInputs = len(givenInputs)
    #inputBinaries = makeInputBin(numInputs)
    initPeriod = 500.0*2**(numInputs)
    period = initPeriod
    if useDefaultInput:
        inputBinaries = makeInputBin(numInputs)
    for i in range(numInputs):
        standardGateName = "IN"+str(i+1)
        #Create the input gate and set the equation to use and its expected protein
        #and promoter truth values which should be the same for an input.
        
        #Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
        tempInput = Gate.Gate(standardGateName,'Input',0.05,2,0.1,13,None)
        if useDefaultInput:
            tempTV = ''.join(inputBinaries[i])
            tempInput.setInputType(['inputs.squInput',period,10,period/2.0,0])
        else:
            tempTV = raw_input(givenInputs[i]+"'s input truthvalue: ")
            tempInput.setInputType(['inputs.specInput',10.0,0.0,initPeriod,tempTV])
        tempInput.setExpectedProtein(list(tempTV))      

        #So the period of the next input alternates twice as frequently as this
        period /= 2.0
        #Add it to the list of inputs. Note the name of the fanOut wire as being
        #from that gate and replace all occurences of that letter with the wire
        #name associated with that gate.
        Inputs.append(tempInput)
        fanOutWireToGateDict[givenInputs[i]] = tempInput
    
    #Makes the intermediate gates
    gateCount = 1
    for i in range(len(givenIntermediates)):
        standardGateName = "G" + str(gateCount)
        givenWireName = givenIntermediates[i]
        if allGatesParsed[givenWireName][0]=="NOR" or allGatesParsed[givenWireName][0]=="NOT": 
            #Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
            tempIntermediate = Gate.Gate(standardGateName,'Repressor',0.05,2,0.1,13,5)
        elif allGatesParsed[givenWireName][0]=="OR" or allGatesParsed[givenWireName][0]=="BUF": 
            tempIntermediate = Gate.Gate(standardGateName,'Activator',0.05,2,0.1,13,5)
        Intermediates.append(tempIntermediate)
        fanOutWireToGateDict[givenWireName] = tempIntermediate
        gateCount += 1
    
    #Makes outputs
    numOutputs = 1
    for i in range(len(givenOutputs)):
        givenWireName = givenOutputs[i]
        standardOutName = "Y"+str(numOutputs)
        if givenWireName==standardOutName:
            standardOutName+="r"
        if allGatesParsed[givenWireName][0]=="NOR" or allGatesParsed[givenWireName][0]=="NOT": 
            standardGateName = "G" + str(gateCount)
            #Gate(name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None)
            tempIntermediate = Gate.Gate(standardGateName,'Repressor',0.05,2,0.1,13,5)
            Intermediates.append(tempIntermediate)
            fanOutWireToGateDict[givenWireName] = tempIntermediate
            allGatesParsed[standardOutName] = ["BUF",givenWireName]      
            
            tempOutput = Gate.Gate(standardOutName,'Output',None,None,None,None,5)
            Outputs.append(tempOutput)
            fanOutWireToGateDict[standardOutName] = tempOutput
            gateCount += 1

        elif allGatesParsed[givenWireName][0]=="OR" or allGatesParsed[givenWireName][0]=="BUF": 
            tempOutput = Gate.Gate(standardOutName,'Output',None,None,None,None,5)
            Outputs.append(tempOutput)
            fanOutWireToGateDict[givenWireName] = tempOutput
        else:
            print "Something isn't right."
        numOutputs += 1

    #Make and connect wires
    wireCount = 1
    wireAliases = {}
    for givenWireName in allGatesParsed.keys():
        #Gets the gate associated with that wire
        tempGate = fanOutWireToGateDict[givenWireName]
        #Ignores the first element which tells the type of gate
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

    return graph, allGates, Inputs, Intermediates, Outputs


#netlist --> graph --> JSON --> General.py opens JSON --> 
def wrapperForNetlist(fileLoc, placeToSave, makeBarGraph=True,makeOtherGraphs=True,useDefaultInput=True,genesToUse=None,Libraries=Libraries):
    """
    Takes in a netlist of a circuit, and a directory to store the intermediate json file.
    Returns the score, if possible, for the arrangement of gates.
    """
    #Make Graph from string (or filepath... checks if fileLoc is a list of strings or single string)
    graph, allGates, Inputs, Intermediates, Outputs = makeGraphFromNetlist(fileLoc,useDefaultInput)
    
    #If we are given gene names to use
    if genesToUse!=None:
        #Throw an error if we are not given the correct number of genes
        if len(genesToUse[0])!=len(Inputs):
            raise IndexError, "Wrong number of input genes. Please try again."
        elif len(genesToUse[1])!=len(Intermediates):
            raise IndexError, "Wrong number of intermediate genes. Please try again."
        elif len(genesToUse[2])!=len(Outputs):
            raise IndexError, "Wrong number of output genes. Please try again."
        
        #Make the lists containing all the geners from the libraries
        allInputs = Optimize.makeInputs(Libraries[0])
        allRepressors = Optimize.makeRepressors(Libraries[1])
        allOutputs = Optimize.makeOutputs(Libraries[2])
        inputsToUse = []
        repressorsToUse = []
        outputsToUse=[]
        #find the genes that were requested
        for geneName in genesToUse[0]:
            for gene in allInputs:
                if gene.getName()==geneName:
                    inputsToUse.append(gene)
                    break
        for geneName in genesToUse[1]:
            for gene in allRepressors:
                if gene.getName()==geneName:
                    repressorsToUse.append(gene)
                    break
        for geneName in genesToUse[2]:
            for gene in allOutputs:
                if gene.getName()==geneName:
                    outputsToUse.append(gene)
                    break
        #The numbers will no longer be the same if one of the genes did not exist in the library
        if len(inputsToUse)!=len(Inputs):
            raise NameError, "One of the inputs did not exist in the library. Please try again."
        elif len(repressorsToUse)!=len(Intermediates):
            raise NameError, "One of the repressors did not exist in the library. Please try again."
        elif len(outputsToUse)!=len(Outputs):
            raise NameError, "One of the outputs did not exist in the library. Please try again."
        
        #swap in the requested genes
        for i in range(len(Inputs)):
            graph.swapGates(Inputs[i],inputsToUse[i])
        for i in range(len(Intermediates)):
            graph.swapGates(Intermediates[i],repressorsToUse[i])
        for i in range(len(Outputs)):
            graph.swapGates(Outputs[i],outputsToUse[i])
        Inputs = inputsToUse
        Intermediates = repressorsToUse
        Outputs = outputsToUse
        
        #check to make sure that no gene was repeated
        if graph.hasRepeatedRepressor():
            raise TypeError, "You cannot use the same gene twice. Please try again."
   
    #don't want to score the circuit if it's sequential
    isSequential = False
    print graph
    if graph.hasLoop():
        isSequential = True
        #So things do not start at equilibrium.
        Intermediates[0].setREUi(1.0)
    
    #Write to Json File (the middleman)
    graph.writeToJson(placeToSave)
    
    #Make graphs from JsonFile
    scores = General.generateDynamicCircuitGraphs(placeToSave, makeBarGraph, makeOtherGraphs, isSequential)
    return scores

        
def makeInputBin(numInputs):
    """
    makes the binary values that will be the default truth values for the inputs
    if the user does not want to specify a wave form.
    """
    
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
    allInputBin.reverse()
    return allInputBin

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