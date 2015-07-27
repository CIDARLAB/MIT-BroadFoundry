# -*- coding: utf-8 -*-
"""
This file is just a scratch paper file for testing out changes to a function 
before putting them in the actual file.
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

    
def makeDAGFromNetlist(netListFileLoc):
    """
    Takes in a circuit in the string form like '((a.b).c)' and returns a dag 
    created from that string, a list of all the gates in that dag, lists of the
    inputs, repressors, and outputs individually from all the gates.
    This assumes all inputs are activators and allows for intermediate genes to
    be activators. Treats all outputs are buffers or ORs.
    """
    
    myFile = open(netListFileLoc,'r')
    netlist = json.load(myFile)
    
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
    allWires = []
    
    for gateInfo in netlist:
        gateInfo = gateInfo.replace(")","")
        gateInfo = gateInfo.replace("(",",")
        temp = gateInfo.split(",")
        allGatesParsed[temp[1]] = []
        allFanOutNames.append(temp[1])
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
    for i in allFanOutNames:
        if i not in allFanInNames:
            givenOutputs.append(i)
    for i in allWireNames:
        print i
        if (i not in givenInputs) and (i not in givenOutputs):
            givenIntermediates.append(i)

    #Make the inputs
    numInputs = len(givenInputs)
    inputBinaries = makeInputBin(numInputs)
    initPeriod = 500.0*2**(numInputs)
    period = initPeriod
    for i in range(numInputs):
        standardGateName = "IN"+str(i+1)
        #Create the input gate and set the equation to use and its expected protein
        #and promoter truth values which should be the same for an input.
        #Gate(name,Km,n,gateType,mB,pB,m_halfLife, p_halfLife)
        tempInput = Gate.Gate(standardGateName,1,2,'Input',30,None,None,None)
        tempInput.setInputType(['inputs.squInput',period,10,period/2.0,0])
        tempInput.setExpectedProtein(inputBinaries[i])
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
    for i in range(len(givenOutputs)):
        givenWireName = givenOutputs[i]
        standardOutName = "Y"+str(numOutputs)
        if allGatesParsed[givenWireName][0]=="NOR" or allGatesParsed[givenWireName][0]=="NOT": 
            standardGateName = "G" + str(gateCount)
            tempIntermediate = Gate.Gate(standardGateName,1000,2,'Repressor',30,20,2,10)
            Intermediates.append(tempIntermediate)
            fanOutWireToGateDict[givenWireName] = tempIntermediate
            allGatesParsed[standardOutName] = ["BUF",givenWireName]      
            
            tempOutput = Gate.Gate(standardOutName,None,None,'Output',None,20,2,10)
            Outputs.append(tempOutput)
            fanOutWireToGateDict[standardOutName] = tempOutput
            gateCount += 1

        elif allGatesParsed[givenWireName][0]=="OR" or allGatesParsed[givenWireName][0]=="BUF": 
            tempOutput = Gate.Gate(standardOutName,None,None,'Output',None,20,2,10)
            Outputs.append(tempOutput)
            fanOutWireToGateDict[givenWireName] = tempOutput
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
    
    dag = DAG.DAG(0,initPeriod,10*initPeriod)
    
    for wire in allWires:
        dag.addWire(wire)
    for gate in allGates:
        dag.addGate(gate)

    return dag, allGates, Inputs, Intermediates, Outputs
    

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



def wrapper(fileLoc):
    """
    Takes in a string representation of a circuit, a set of libraries of inputs,
    repressors, and outputs, and a directory to store the intermediate json file.
    Prints the graphs generated for mRNA and protein concentrations generated
    by General.
    """
    #Make DAG from string
    dag, allGates, Inputs, Intermediates, Outputs = makeDAGFromNetlist(fileLoc)


    print dag
    #Write to Json File
    #fileLoc = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Example5.json'
#    dag.writeToJson(fileLoc)
    #Make graphs from JsonFile
#    General.generateDynamicCircuitGraphs(fileLoc, makeBarGraph)
    
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