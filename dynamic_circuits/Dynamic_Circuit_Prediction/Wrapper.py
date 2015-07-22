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
    if inputNames[0] == "0":
        inputNames2 = inputNames[1:]
    elif inputNames[0] != "0":
        inputNames2 = inputNames[:]
    
    inputBin = makeInputBin(len(inputNames2))
    
    #Make a gate for all but 0
    initPeriod = 500*2**(len(inputNames2))
    period = initPeriod
    for i in range(len(inputNames2)):
        name = inputNames2[i]
        tempInput = Gate.Gate(name,10,2,'Input',30,None,None,None)
        tempInput.setInputType(['inputs.squInput',period,100,0,1])
        tempInput.setExpectedProtein(inputBin[i])
        period /= 2.0 
        Inputs.append(tempInput)
#    print stringInputs
#    print inputNames
#    print inputNames2
#    pause = raw_input("pause1")
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
#    print circuitString
#    print "Input fan out"
#    print inputFanOutDict
#    pause = raw_input("pause2")
    similarWires = {}
    for i in inputFanOutDict.keys():
        sameWires = inputFanOutDict[i]
        if len(sameWires)>1:
            mainWire = sameWires[0]
            similarWires[mainWire] = sameWires
            for wire in sameWires[1:]:
                circuitString = circuitString.replace(wire,mainWire)
    
#    print similarWires
#    print circuitString
#    pause = raw_input("pause3")
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
        fanIn = piece[1:-1].split(".")
        for wire in inputFanOutDict["0"]:
            try:
                fanIn.remove(wire)
            except ValueError:
                pass
        allGatesDictFanIn["G"+str(gateCount)] = fanIn
        
#        print circuitString
#        print allGatesFanOut
#        print allGatesDictFanIn
#        print gateCount
        
        gateCount += 1
        
#        pause = raw_input("pause4")
        
    outputDict["Y"] = [circuitString]
#    print outputDict
#    print circuitString
#    print allGatesFanOut
#    print allGatesDictFanIn
#    print gateCount
#    pause = raw_input("pause5")
        
#    print "Repressor fan in"    
#    print allGatesDictFanIn
#    print "Repressor fan out"
#    print allGatesFanOut
#    print "Output fan in"
#    print outputDict
#    print "Wires from the same place"
#    print similarWires
#    pause = raw_input("pause6")
    allWireNames = []
    whichGaveWhich = []
    for gate in allGatesDictFanIn:
        for wire in allGatesDictFanIn[gate]:
            allWireNames.append(wire)
            whichGaveWhich.append(gate)
    for otpt in outputDict:
        for wire in outputDict[otpt]:
            allWireNames.append(wire)
            whichGaveWhich.append(otpt)
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
            while newWireName in aliases:
                newWireName = groupOfSim[groupOfSim.index(newWireName)+1]
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
        tempRepressor = Gate.Gate(repressorName,1000,2,'Repressor',30,20,2,10)
        Repressors.append(tempRepressor)
    allOutputNames = outputDict.keys()
    Outputs = []
    for outputName in allOutputNames:
        tempOut = Gate.Gate(outputName,None,None,'Output',None,20,2,10)
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
    dag = DAG.DAG(0,initPeriod,10*initPeriod)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
#    pause = raw_input("pause7")
#    print dag
    return dag, allGates, Inputs, Repressors, Outputs
    
    
def makeGatesFromLibraries(Libraries,allGates,Inputs,Repressors,Outputs):
    inputsDir = Libraries[0]
    repressorsDir = Libraries[1]
    outputsDir = Libraries[2]
    myFile = open(inputsDir,'r')
    inputsFromFile = json.load(myFile)
    myFile.close()
    myFile = open(repressorsDir,'r')
    repressorsFromFile = json.load(myFile)
    myFile.close()
    myFile = open(outputsDir,'r')
    outputsFromFile = json.load(myFile)
    myFile.close()
    
    #Make Inputsr,Outputsr, Repressorsr, allGatesr
    Inputsr = []
    Repressorsr = []
    Outputsr = []
    for i in range(len(Inputs)):
        gateInfo = inputsFromFile[i]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),gateInfo["Km"],gateInfo["n"],"Input",gateInfo["mB"],None,None,None)
        Inputsr.append(tempGate)
    for i in range(len(Repressors)):
        gateInfo = repressorsFromFile[i]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),gateInfo["Km"],gateInfo["n"],"Repressor",gateInfo["mB"],gateInfo["pB"],gateInfo["m_halfLife"],gateInfo["p_halfLife"])
        Repressorsr.append(tempGate)
    for i in range(len(Outputs)):
        gateInfo = outputsFromFile[i]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),None,None,"Output",None,gateInfo["pB"],gateInfo["m_halfLife"],gateInfo["p_halfLife"])
        Outputsr.append(tempGate)

    allGatesr = Inputsr + Repressorsr + Outputsr
    return allGatesr, Inputsr, Repressorsr, Outputsr
    
def performSwaps(dag,allGates,Inputs,allGatesr,Inputsr):
    #Make swaps in dag    
    for i in range(len(allGates)):       
        if allGatesr[i].getGateType() != "Output":
            allGatesr[i].setFanOut(allGates[i].getFanOut())
        if allGatesr[i].getGateType() != "Input":
            allGatesr[i].setFanIn(allGates[i].getFanIn())
    for i in range(len(Inputsr)):
        Inputsr[i].setInputType(Inputs[i].getInputType())
        Inputsr[i].setExpectedProtein(Inputs[i].getExpectedProtein())
    #add new gates to dag 
    for gate in allGatesr:
        dag.addGate(gate)
    #remove old gates from dag
    for gate in allGates:
        dag.removeGate(gate)

circuitString1 = '((((a.0).(b.0)).0).(c.0))' #00000001
circuitString2 = '((((a.b).(a.c)).a).(b.c))' #00010111
circuitString3 = '(((((a.c).c).0).(b.c)).(((a.c).c).b))' #00111001
circuitString4 = 'a' #00001111
circuitString5 = '(a.b)' #11000000
circuitString6 = '((a.b).0)' #00111111
circuitString7 = '(a.0)' #11110000
circuitString8 = '((a.b).c)' #00101010
circuitString9 = '((a.0).a)' #00000000
circuitString10 = '((a.0).b)' #00001100

fileLoc = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015/test.json"
inputsDir = "C:\Users\Arinze\Documents\GitHub\MIT-BroadFoundry\dynamic_circuits\DynamicCircuits\Libraries/InputLibrary2.json"
repressorsDir = "C:\Users\Arinze\Documents\GitHub\MIT-BroadFoundry\dynamic_circuits\DynamicCircuits\Libraries/RepressorLibrary2.json"
outputsDir = "C:\Users\Arinze\Documents\GitHub\MIT-BroadFoundry\dynamic_circuits\DynamicCircuits\Libraries/OutputLibrary2.json"
Libraries = [inputsDir, repressorsDir, outputsDir]
def wrapper(circuitString, Libraries, fileLoc):
    #Make DAG from string
    dag, allGates, Inputs, Repressors, Outputs = makeDAGFromString(circuitString)
    #Get replacements for the gates
    allGatesr, Inputsr, Repressorsr, Outputsr = makeGatesFromLibraries(Libraries,allGates, Inputs, Repressors, Outputs)
    #Perform Swaps
#    performSwaps(dag, allGates, Inputs, allGatesr, Inputsr)

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

def func(circuitString):
    """
    pausing in between each set of libraries, generates a graph for each 
    library for a given circuit.
    """
    for i in range(1,10):
         inputsDir = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Libraries/InputLibrary"+str(i)+".json"
         repressorsDir = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Libraries/RepressorLibrary"+str(i)+".json"
         outputsDir = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Libraries/OutputLibrary"+str(i)+".json"
         Libraries = [inputsDir, repressorsDir, outputsDir]
         wrapper(circuitString, Libraries, fileLoc)
         pause = raw_input("Finished "+str(i))
         
def makeInputBin(numInputs):
    binSize = 2**numInputs
    base = "10"
    allInputBin = []
    for i in range(numInputs):
        
        tempInput = base
        while len(tempInput)<binSize:
            tempInput += base
        nextLen = len(base)*2
        while len(base)<nextLen:
            base = "1" + base + "0"
        allInputBin.append(list(tempInput))
        #allInputBin.append(tempInput)
    allInputBin.reverse()
    return allInputBin