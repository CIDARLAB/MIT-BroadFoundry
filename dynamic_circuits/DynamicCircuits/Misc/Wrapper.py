# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 11:45:31 2015

@author: Arinze
"""

import DAG
import Gate
import Wire
import inputs
import json
import General
import GeneralJsonIO
import inputs
import update
import re
circuitString = "(((((a.0).b).(a.0)).0).(c.0))"
circuitString2 = "(((((a.0).b).(a.0)).0).(c.0))"
circuitString3 = "(((((((b.0).a).(b.0)).c).(((b.0).a).(b.0))).((b.0).c)).0)"
def makeDAGFromString(circuitString):
    #Given a string representation of a circuit and the Input pattern.
    #Generate the DAG with 0 values for gate properties.
    
    #Splt the string up into the inputs and then remove duplicates.
    print circuitString
    allInputs = re.findall(r"[\w]", circuitString)
    stringInputs = []
    inputNames = []
    inputCount = 1
    inputFanOutDict = {}
    #Replace the input name with IN then a number
    
    for i in allInputs:
        if i not in stringInputs:
            stringInputs.append(i)
            inputNames.append("IN"+str(inputCount))
            inputCount += 1
    stringInputs.sort()


    #Create the list of Input gates
    Inputs = []
    if stringInputs[0]=="0":
        inputNames = inputNames[:-1]
        inputNames = ["0"]+inputNames
    #Make a gate for all but 0
    for name in inputNames[1:]:
        tempInput = Gate.Gate(name,0,0,'Input',0,None,None,None)
        Inputs.append(tempInput)
    print stringInputs
    print inputNames
    wireCount = 1
    inputFanOutDict["0"] = []
    for i in range(len(inputNames)):
#        print circuitString
        inputFanOutDict[inputNames[i]] = []
        while circuitString.count(stringInputs[i])!=0:
            circuitString = circuitString.replace(stringInputs[i],"W"+str(wireCount),1)
            inputFanOutDict[inputNames[i]].append("W"+str(wireCount))
            wireCount += 1
#            print circuitString
    print "Input fan out"
    print inputFanOutDict
    similarWires = {}
    for i in inputFanOutDict.keys():
        sameWires = inputFanOutDict[i]
        mainWire = sameWires[0]
        similarWires[mainWire] = sameWires
        for wire in sameWires[1:]:
            circuitString = circuitString.replace(wire,mainWire)
    parIndex = findSmallestParentheses(circuitString)
    allGatesDictFanIn = {}
    allGatesFanOut = {}
    outputDict = {}
    gateCount = 1
    while parIndex != None:
        
        piece = circuitString[parIndex[0]:parIndex[1]+1]
        allGatesFanOut["G"+str(gateCount)] = []
        while circuitString.count(piece)!=0:
            circuitString = circuitString.replace(piece,"W"+str(wireCount),1)
            allGatesFanOut["G"+str(gateCount)].append("W"+str(wireCount))
            wireCount += 1
            
        sameWires = allGatesFanOut["G"+str(gateCount)]
        mainWire = sameWires[0]
        similarWires[mainWire] = sameWires
        for wire in sameWires[1:]:
            circuitString = circuitString.replace(wire,mainWire)
            
        
        parIndex = findSmallestParentheses(circuitString)
        #A zero as a fan in is the same as a not/a single fan in.
        if parIndex != None:
            fanIn = piece[1:-1].split(".")
            for wire in inputFanOutDict["0"]:
                try:
                    fanIn.remove(wire)
                except ValueError:
                    pass
            allGatesDictFanIn["G"+str(gateCount)] = fanIn
        else:
            allGatesFanOut.pop("G"+str(gateCount))
            fanIn = piece[1:-1].split(".")
            for wire in inputFanOutDict["0"]:
                try:
                    fanIn.remove(wire)
                except ValueError:
                    pass
            outputDict["Y"] = fanIn
        
        gateCount += 1
    print "Repressot fan in"    
    print allGatesDictFanIn
    print "Repressor fan out"
    print allGatesFanOut
    print "Output fan in"
    print outputDict
    print "Wires from the same place"
    print similarWires
    allWireNames = []
    whichGaveWhich = []
    for gate in allGatesDictFanIn:
        for wire in allGatesDictFanIn[gate]:
            allWireNames.append(wire)
            whichGaveWhich.append(gate)
    for otpt in outputDict:
        for wire in outputDict[otpt]:
            allWireNames.append(wire)
            whichGaveWhich.append(gate)
    mapToMainWire = {}
    for mainWire in similarWires:
        for wire in similarWires[mainWire]:
            mapToMainWire[wire] = mainWire

    #Rename wires and in gates' fanIn and fanOut to prevent the same wire having
    #two destinations. Use alias to rename wires so wire numbers don't skip
    #Remove ghost wires in the process.
    wireCount = 1
    allWireNames2 = allWireNames[:]
    allWireNames = []
    aliases = {}
    for i in range(len(allWireNames2)):
        if allWireNames2[i] not in aliases:
            allWireNames.append("W"+str(wireCount))
            aliases[allWireNames2[i]] = "W"+str(wireCount)
        elif allWireNames2[i] in aliases:
            oldWireName = allWireNames2[i]
            mainWire = mapToMainWire[oldWireName]
            groupOfSim = similarWires[mainWire]
            newWireName = groupOfSim[groupOfSim.index(oldWireName)+1]
            allWireNames.append("W"+str(wireCount))
            aliases[newWireName] = "W"+str(wireCount)
            print oldWireName,"changed to", newWireName
            try:
                allGatesDictFanIn[whichGaveWhich[i]].remove(oldWireName)
                allGatesDictFanIn[whichGaveWhich[i]].append(newWireName)
            except KeyError:
                outputDict[whichGaveWhich[i]].remove(oldWireName)
                outputDict[whichGaveWhich[i]].append(newWireName)
        wireCount += 1
    
    inputFanOutDict.pop("0")
    allDicts = [inputFanOutDict,allGatesDictFanIn,allGatesFanOut,outputDict]
    allDicts2 = []
    for fanDict in allDicts:
        tempDict = {}
        for key in fanDict.keys():
            gateFan = fanDict[key]
            tempDict[key] = []
            for wire in gateFan:
                try:
                    tempDict[key].append(aliases[wire])
                #Skip ghost wires
                except KeyError:
                    pass
        allDicts2.append(tempDict)
    inputFanOutDict = allDicts2[0]
    allGatesDictFanIn = allDicts2[1]
    allGatesFanOut = allDicts2[2]
    outputDict = allDicts2[3]
    
    allWires = []
    for wireName in allWireNames:
        tempWire = Wire.Wire(wireName)
        allWires.append(tempWire)
    allRepressorNames = allGatesDictFanIn.keys()
    Repressors = []
    for repressorName in allRepressorNames:
        tempRepressor = Gate.Gate(repressorName,0,0,'Repressor',0,0,0,0)
        Repressors.append(tempRepressor)
    allOutputNames = outputDict.keys()
    Outputs = []
    for outputName in allOutputNames:
        tempOut = Gate.Gate(outputName,None,None,'Output',None,0,0,0)
        Outputs.append(tempOut)
    allGates = Inputs + Repressors + Outputs

    #Attach wires utilizing the fact that allWires is in same order as allWireNames
    
    #Make lists that are the same order as Inputs, Repressors and Outputs but
    #have the names for easy indexing.
    inputNames = []
    for ipt in Inputs:
        inputNames.append(ipt.getName())
        
    repressorNames = []
    for rep in Repressors:
        repressorNames.append(rep.getName())
    outputNames = []
    for otpt in Outputs:
        outputNames.append(otpt.getName())
    #Attach inputs to their fanOut
    for i in range(len(inputNames)):
        ipt = Inputs[i]
        for wireName in inputFanOutDict[inputNames[i]]:
            wireIndex = allWireNames.index(wireName)
            ipt.addFanOutWire(allWires[wireIndex])
    #Attach repressors to their fanIn
    for i in range(len(repressorNames)):
        rep = Repressors[i]
        for wireName in allGatesDictFanIn[repressorNames[i]]:
            wireIndex = allWireNames.index(wireName)
            rep.addFanInWire(allWires[wireIndex])
    #Attach repressors to their fanOut
    for i in range(len(repressorNames)):
        rep = Repressors[i]
        for wireName in allGatesFanOut[repressorNames[i]]:
            wireIndex = allWireNames.index(wireName)
            rep.addFanOutWire(allWires[wireIndex])
    #Attach outputs to their fanIn
    for i in range(len(outputNames)):
        otpt = Outputs[i]
        for wireName in outputDict[outputNames[i]]:
            wireIndex = allWireNames.index(wireName)
            otpt.addFanInWire(allWires[wireIndex])
    #Set input types 
#    a.setInputType(['inputs.squInput',1000,100,0,0])    
#    b.setInputType(['inputs.squInput',500,100,0,0])
    dag = DAG.DAG(0,2000,500)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
    print dag
    return dag, allGates, Inputs, Repressors, Outputs
    
    
def makeGatesFromLibraries(Libraries,allGates,Inputs,Repressors,Outputs):
    #Make Inputsr,Outputsr, Repressorsr, allGatesr
    ar = Gate.Gate('a',40,2,'Input',30,None,None,None)
    br = Gate.Gate('b',40,2,'Input',30,None,None,None)
    G1r = Gate.Gate('G1',30000,2,'Repressor',30,20,2,10)
    G2r = Gate.Gate('G2',30000,2,'Repressor',30,20,2,10)
    yr = Gate.Gate('y',None,None,'Output',None,20,2,10)
    Inputsr = [ar,br]
    Repressorsr = [G1r,G2r]
    Outputsr = [yr]
    allGatesr = Inputsr + Repressorsr + Outputsr
    return allGatesr, Inputsr, Repressorsr, Outputsr
    
def performSwaps(dag,allGates,Inputs,Repressors,Outputs,allGatesr,Inputsr,RepressorsR,Outputsr):
    #Make swaps in dag    
    for i in range(len(allGates)):
        allGatesr[i].setFanOut(allGates[i].getFanOut())
        allGatesr[i].setFanIn(allGates[i].getFanIn())
    for i in range(len(Inputsr)):
        Inputsr[i].setInputType(Inputs[i].getInputType())
    #add new gates to dag 
    for gate in allGatesr:
        dag.addGate(gate)
    #remove old gates from dag
    for gate in allGates:
        dag.removeGate(gate)

def wrapper(circuitString, Libraries, fileLoc):
    #Make DAG from string
    dag, allGates, Inputs, Repressors, Outputs = makeDAGFromString(circuitString)
    #Get replacements for the gates
    allGatesr, Inputsr, Repressorsr, Outputsr = makeGatesFromLibraries(Libraries,allGates, Inputs, Repressors, Outputs)
    #Perform Swaps
    performSwaps(dag, allGates, Inputs, Repressors, Outputs, allGatesr, Inputsr, Repressorsr, Outputsr)

    print dag
    #Write to Json File
    #fileLoc = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Example5.json'
    dag.writeToJson(fileLoc)
    #Make graphs from JsonFile
    General.generateDynamicCircuitGraphs(fileLoc)
    
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