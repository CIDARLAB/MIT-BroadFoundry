# -*- coding: utf-8 -*-


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
import matplotlib.pyplot as plt
import numpy as np

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

    if inputNames[0]=="0":
        inputNames2 = inputNames[1:]
    elif inputNames[0]!="0":
        inputNames2 = inputNames[:]
    initPeriod = 5000*2**(len(inputNames2))
    period = 5000*2**(len(inputNames2))
    #Make a gate for all but 0
    for name in inputNames2:
        tempInput = Gate.Gate(name,40,2,'Input',30,None,None,None)
        tempInput.setInputType(['inputs.squInput',period,100,0,0])
        period /= 2.0
        Inputs.append(tempInput)
    wireCount = 1
    if inputNames[0] == "0":
        inputFanOutDict["0"] = []
    for i in range(len(inputNames)):
#        print circuitString
        inputFanOutDict[inputNames[i]] = []
        while circuitString.count(stringInputs[i])!=0:
            circuitString = circuitString.replace(stringInputs[i],"W"+str(wireCount),1)
            inputFanOutDict[inputNames[i]].append("W"+str(wireCount))
            wireCount += 1
#            print circuitString
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
            try:
                for wire in inputFanOutDict["0"]:
                    try:
                        fanIn.remove(wire)
                    except ValueError:
                        pass
            except KeyError:
                pass
            allGatesDictFanIn["G"+str(gateCount)] = fanIn
        else:
            allGatesFanOut.pop("G"+str(gateCount))
            fanIn = piece[1:-1].split(".")
            try:
                for wire in inputFanOutDict["0"]:
                    try:
                        fanIn.remove(wire)
                    except ValueError:
                        pass
            except KeyError:
                pass
            outputDict["Y"] = fanIn
        
        gateCount += 1
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
    try:
        inputFanOutDict.pop("0")
    except KeyError:
        pass
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
        tempRepressor = Gate.Gate(repressorName,30000,2,'Repressor',30,20,2,10)
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
    dag = DAG.DAG(0,initPeriod,10000)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
    return dag, allGates, Inputs, Repressors, Outputs
    
def makeDAGFromString2(circuitString):
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
    

    
    #Make a gate for all but 0
    initPeriod = 500*2**(len(inputNames2))
    period = initPeriod
    for name in inputNames2:
        tempInput = Gate.Gate(name,40,2,'Input',30,None,None,None)
        tempInput.setInputType(['inputs.squInput',period,100,0,0])
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
        tempRepressor = Gate.Gate(repressorName,30000,2,'Repressor',30,20,2,10)
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
    print dag
    return dag, allGates, Inputs, Repressors, Outputs
    
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
fileLoc1 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test1.json"
fileLoc2 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test2.json"
fileLoc3 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test3.json"
fileLoc4 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test4.json"
fileLoc5 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test5.json"
fileLoc6 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test6.json"
fileLoc7 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test7.json"
fileLoc8 = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Test/test8.json"

def wrapper(circuitString, fileLoc=""):
    #Make DAG from string
#    dag, allGates, Inputs, Repressors, Outputs = makeDAGFromString(circuitString)
#    print dag
    
    dag2, allGates2, Inputs2, Repressors2, Outputs2 = makeDAGFromString2(circuitString)
#    print dag2
    #Write to Json File
    #fileLoc = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Example5.json'
#    dag.writeToJson(fileLoc)
    dag2.writeToJson(fileLoc)
#    #Make graphs from JsonFile
    General.generateDynamicCircuitGraphs(fileLoc)
    
values = [0.5,5,4,0.8]
tv = [0,1,1,0] 
   
def func3(values,tv):
#    tzero = [1,4]
#    yzero = [0.5,0.8]
#    tone = [2,3]
#    yone = [5,4]
#    
#    zeros = (tzero,yzero)
#    ones = (tone,yone)
#        # the width of the bars: can also be len(x) sequence
#
#    
#    combinedOne = [0,5,4,0]
#    combinedZero = [0.5,0,0,0.8]
#    width = 0.425  
#    t = np.arange(4)
#    plt.figure()
#    #Plot each protein against time with its name as the label
#    plt.axis(xmin=0, xmax=6, ymin=0,ymax=6)
#    
#    plt.bar(t,combinedOne,color="b",label = "ones")
#    plt.bar(t,combinedZero,color="r",label = "zeros")
#    plt.xticks(t+width, ('11', '10', '01', '00') )  
#    plt.xlabel('Time (min)')
#    plt.ylabel('Proteins per cell')
#    plt.title('Concentration of Protein Over Time')
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    plt.show()
    numInputs = 2
    numBins = 2**numInputs
    binVals = getBinaryInReverse(numBins)
    width = 0.425  
    t = np.arange(numBins)
    plt.figure()
    #Plot each protein against time with its name as the label
    plt.axis(xmin=0, xmax=6, ymin=0,ymax=6)
    barlist = plt.bar(range(numBins),values)
    for i in range(numBins):
        if tv[i]==0:
            barlist[i].set_color('r')
    plt.xticks(t+width, binVals )  
    plt.xlabel('Time (min)')
    plt.ylabel('Proteins per cell')
    plt.title('Concentration of Protein Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
def getBinaryInReverse(maxVal):
    binVals = []
    decVals = range(maxVal)
    decVals.reverse()
    for i in decVals:
        temp="{0:b}".format(i)
        if len(binVals) == 0:
            binVals.append(temp)
        else:
            while len(temp)<len(binVals[0]):
                temp = "0"+temp
            binVals.append(temp)
    return binVals