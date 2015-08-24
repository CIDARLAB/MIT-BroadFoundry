# -*- coding: utf-8 -*-
"""

Used for creating graphs from netlists and checking if they are sequential and
if they stabilize to any point.

"""

import json
import Objects
#import LogicOperations

#Latches
Gated_D_Latch_FileLoc = "JsonFiles/SR_Latches/Gated_D_Latch.json"
Gated_D_Latch_placeToSave = "JsonFiles/SR_Latches/Gated_D_Latch_DAG.json"
Gated_D_Latch_noCInv_FileLoc = "JsonFiles/SR_Latches/Gated_D_Latch_noCInv.json"
Gated_D_Latch_noCInv_placeToSave = "JsonFiles/SR_Latches/Gated_D_Latch_noCInv_DAG.json"
Negative_Edge_Triggered_D_Flip_Flop_FileLoc = "JsonFiles/SR_Latches/Negative_Edge_Triggered_D_Flip_Flop.json"
Negative_Edge_Triggered_D_Flip_Flop_placeToSave = "JsonFiles/SR_Latches/Negative_Edge_Triggered_D_Flip_Flop_DAG.json"
Transparent_D_Latch_noDinv_FileLoc = "JsonFiles/SR_Latches/Transparent_D_Latch_noDinv.json"
Transparent_D_Latch_noDinv_placeToSave = "JsonFiles/SR_Latches/Transparent_D_Latch_noDinv_DAG.json"
Transparent_D_Latch_noInpinv_FileLoc = "JsonFiles/SR_Latches/Transparent_D_Latch_noInpinv.json"
Transparent_D_Latch_noInpinv_placeToSave = "JsonFiles/SR_Latches/Transparent_D_Latch_noInpinv_DAG.json"
#TestCase
Test1_FileLoc = "JsonFiles/SR_Latches/Test1.json"
Test1_placeToSave = "JsonFiles/SR_Latches/Test1_DAG.json"
Test2_FileLoc = "JsonFiles/SR_Latches/Test2.json"
Test2_placeToSave = "JsonFiles/SR_Latches/Test2_DAG.json"

def makeGraphFromNetlist(netlist):
    if type(netlist)==str:
        myFile = open(netlist,'r')
        netlist = json.load(myFile)
        myFile.close()
    #Initialize the lists of components
    Inputs = []
    Gates = []
    
    givenNameToGateDict = {}
    allGatesParsed = {}
    allFanOutNames = []
    allFanInNames = []
    allGateNames = []
    
    for gateInfo in netlist:
        gateInfo = str(gateInfo)
        gateInfo = gateInfo.replace(")","")
        gateInfo = gateInfo.replace("(",",")
        temp = gateInfo.split(",")
        allGatesParsed[temp[1]] = []
        allFanOutNames.append(temp[1])
        if temp[1] not in allGateNames:
            allGateNames.append(temp[1])
        for i in range(len(temp)):
            if i!=1:
                allGatesParsed[temp[1]].append(temp[i])
            if i>1:
                if temp[i] not in allFanInNames:
                    allFanInNames.append(temp[i])
                if temp[i] not in allGateNames:
                    allGateNames.append(temp[i])
             
    givenInputs = []
    givenGates = []
    for i in allFanInNames:
        if i not in allFanOutNames:
            givenInputs.append(i)
    givenInputs.sort()
    
    for i in allGateNames:
        if i not in givenInputs:
            givenGates.append(i)
            
    #Make the inputs
    numInputs = len(givenInputs)
    inputAliases = {}
    for i in range(numInputs):
        standardGateName = "IN"+str(i+1)
        #Create the input gate
        #Gate(name,gateType,state="0")
        startState = -1
        while startState not in ["1","0"]:
            print "A start state can only be 1 or 0."
            startState = raw_input("Enter the start state for "+givenInputs[i]+": ")
            print ""
        tempInput = Objects.Gate(standardGateName,'INPUT',startState)
        inputAliases[givenInputs[i]]=standardGateName

        #Add it to the list of inputs. Note the name of the fanOut wire as being
        #from that gate and replace all occurences of that letter with the wire
        #name associated with that gate.
        Inputs.append(tempInput)
        givenNameToGateDict[givenInputs[i]] = tempInput
    
    #Make the gates
    gateCount = 1
    gateAliases = {}
    for i in range(len(givenGates)):
        standardGateName = "G" + str(gateCount)
        givenName = givenGates[i]
        startState = -1
        while startState not in ["1","0"]:
            print "A start state can only be 1 or 0."
            startState = raw_input("Enter the start state for "+givenName+": ")
            print ""
        tempGate = Objects.Gate(standardGateName,allGatesParsed[givenName][0],startState)
        Gates.append(tempGate)
        givenNameToGateDict[givenName] = tempGate
        gateAliases[givenName] = standardGateName
        gateCount += 1
    

    #Connect gates
    for givenName in allGatesParsed.keys():
        #get the gate associated with that wire
        tempGate = givenNameToGateDict[givenName]
        #Ignore the first element which tells the type of gate
        fanInNames = allGatesParsed[givenName][1:]
        for gateName in fanInNames:
            tempGate.addFanIn(givenNameToGateDict[gateName])

    allGates = Inputs + Gates
    graph = Objects.Graph()
    for gate in allGates:
        graph.addGate(gate)
        
    print inputAliases
    print gateAliases
    
#    print graph
#    print graph.isSequential()
    return graph, allGates

def test():

    graph = makeGraphFromNetlist(Gated_D_Latch_FileLoc)[0]
    print graph.isSequential()
    print graph
    print "break"
    graph.printHistory()
    print "break"
    graph.runUntilStableOrRepeat(0)

    
def test2():

    graph = makeGraphFromNetlist(Gated_D_Latch_noCInv_FileLoc)[0]
    print graph.isSequential()
    print graph
    print "break"
    graph.runUntilStableOrRepeat(0)
    print "break"
