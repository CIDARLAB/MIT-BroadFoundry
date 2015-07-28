# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 12:09:16 2015

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
import itertools

inputsDir = "JsonFiles/Libraries/InputLibrary2.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary2Edit.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary2.json"
Libraries = [inputsDir, repressorsDir, outputsDir]
    
repressilatorFileLoc = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan/repressilator.json"
placeToSaveRepressilator = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan/repressilatorDAG.json"
#
truthValueExampleFileLoc = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan\SplitByTruthValue_OR/01101001.json" #01101001
placeToSaveTruthValueExample = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan/01101001DAG.json"
netlistLoc = "JsonFiles/standardOR.json"
intermediateFile = "JsonFiles/intermediateDAG.json"

truthValue = "00010001"
inputNames = ["IPTG1","IPTG2","IPTG3"]
outputNames = ["FP1"]

def makeDAGFromNetlist(inputNetlist):
    """
    Takes in a circuit in the string form like '((a.b).c)' and returns a dag 
    created from that string, a list of all the gates in that dag, lists of the
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
    for i in allFanOutNames:
        if i not in allFanInNames:
            givenOutputs.append(i)
    for i in allWireNames:
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
    
    if numInputs > 0:
        dag = DAG.DAG(0,initPeriod,initPeriod)
    else:
        #arbitrary
        dag = DAG.DAG(0,1000,10000)
    for wire in allWires:
        dag.addWire(wire)
    for gate in allGates:
        dag.addGate(gate)

    return dag, allGates, Inputs, Intermediates, Outputs

    
def performSwaps(dag,allGates,Inputs,allGatesr,Inputsr):
    """
    Given a dag, a set of gates in the dag and a new set of gates, it will swap
    in the new gates for the old gates. The number of each type of gate must be
    the same.
    """
    
    #Make swaps in dag
#    print "Performing swaps"
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

def wrapperForNetlist(fileLoc, placeToSave, makeBarGraph=True):
    """
    Takes in a string representation of a circuit, a set of libraries of inputs,
    repressors, and outputs, and a directory to store the intermediate json file.
    Prints the graphs generated for mRNA and protein concentrations generated
    by General.
    """
    #Make DAG from string
    dag, allGates, Inputs, Intermediates, Outputs = makeDAGFromNetlist(fileLoc)


    print dag
    if dag.hasLoop():
        makeBarGraph = False
        #So things do not start at equilibrium.
        Intermediates[0].setInitialProtein(1.0)
    #Write to Json File
    dag.writeToJson(placeToSave)
    #Make graphs from JsonFile
    General.generateDynamicCircuitGraphs(placeToSave, makeBarGraph)
    
         
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
        dag, allGates, Inputs, Intermediates, Outputs = makeDAGFromNetlist(netlist)
        print dag
        IntermediatesrList = makeReressorsGroups(Libraries[1],len(Intermediates))
        for IntermediatesrTuple in IntermediatesrList:
            #-------------------------------
            dag, allGates, Inputs, Intermediates, Outputs = makeDAGFromNetlist(netlist)
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
            performSwaps(dag, allGates, Inputs, allGatesr, Inputsr)
#            print dag
            dag.writeToJson(intermediateFile)
            if dag.hasLoop():
                Intermediatesr[0].setInitialProtein(1.0)
                makeBarGraph = False
                
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
    return bestNetlist, gateAliases
            

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
    
def makeReressorsGroups(repressorLibrary,numRepressors):
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
