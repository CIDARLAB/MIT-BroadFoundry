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
        dag = DAG.DAG(0,initPeriod,10*initPeriod)
    else:
        #arbitrary
        dag = DAG.DAG(0,1000,10000)
    for wire in allWires:
        dag.addWire(wire)
    for gate in allGates:
        dag.addGate(gate)

    return dag, allGates, Inputs, Intermediates, Outputs

def makeDAGFromString(circuitString):
    """
    Takes in a circuit in the string form like '((a.b).c)' and returns a dag 
    created from that string, a list of all the gates in that dag, lists of the
    inputs, repressors, and outputs individually from all the gates.
    """
    print circuitString
    
    #Initialize the lists of components
    Inputs = []
    Repressors = []
    Outputs = []
    allWires = []
    #wireCount is used to name wires
    wireCount = 1
    #Maps wire names to the gate they come from
    fanOutWireToGateDict = {}
    
    #Remove the zeros because they do not get a wire, but leave the parentheses
    #so we still know to invert.
    circuitString = circuitString.replace(".0","")

    #Find the number of unique inputs and sort them so a is the most significant
    #and c is least significant
    allInputsWithDuplicates = re.findall(r"[\w]", circuitString)
    allInputLetters = []
    for i in allInputsWithDuplicates:
        if i not in allInputLetters:
            allInputLetters.append(i)
    allInputLetters.sort()
    
    #Get the number of inputs
    numInputs = len(allInputLetters)
    #Find the binary values for each input. This is a list with the truth
    #value for the most significant bit first.
    inputBinaries = makeInputBin(numInputs)

    #Make the inputs
    #We want will want the period of each input to be a constant multiplied by
    #a power of two
    initPeriod = 500.0*2**(numInputs)
    period = initPeriod
    for i in range(numInputs):
        #Set the input name and wire name
        gateName = 'IN'+str(i+1)
        wireName = 'W' + str(wireCount)
        #Create the input gate and set the equation to use and its expected protein
        #and promoter truth values which should be the same for an input.
        #Gate(name,Km,n,gateType,mB,pB,m_halfLife, p_halfLife)
        tempInput = Gate.Gate(gateName,1,2,'Input',30,None,None,None)
        tempInput.setInputType(['inputs.squInput',period,10,period/2.0,0])
        tempInput.setExpectedProtein(inputBinaries[i])
        #So the period of the next input alternates twice as frequently as this
        period /= 2.0
        #Add it to the list of inputs. Note the name of the fanOut wire as being
        #from that gate and replace all occurences of that letter with the wire
        #name associated with that gate.
        Inputs.append(tempInput)
        fanOutWireToGateDict[wireName] = tempInput
        circuitString = circuitString.replace(allInputLetters[i],wireName)
        wireCount += 1
    
    #Make the gates
    gateCount = 1
    parIndex = findSmallestParentheses(circuitString)
    while parIndex != None:
        gateName = "G"+str(gateCount)
        wireName = "W" + str(wireCount)
        #Get the string of the subcircuit and use it to determine which wires
        #should be in the fanIn
        subCircuit = circuitString[parIndex[0]:parIndex[1]+1]
        fanIn = subCircuit[1:-1].split(".")
        #Make a the gate, add it to the list of repressors, and associate the 
        #fanOut wire name with the gate.
        #Gate(name,Km,n,gateType,mB,pB,m_halfLife, p_halfLife)
        tempRepressor = Gate.Gate(gateName,1000,2,'Repressor',30,20,2,10)
        Repressors.append(tempRepressor)
        fanOutWireToGateDict[wireName] = tempRepressor
        #Replace all occurrences of subcircuit in the circuitString with wireName.
        circuitString = circuitString.replace(subCircuit,wireName)
        #Create wires for each of the fanIn of the gate we just made. Attach 
        #them to the gate itself and to the gate that the wire comes from using
        #the dictionary mapping wire names to gates to find which gate the wire
        #comes from.
        for fanInWireName in fanIn:
            wire = Wire.Wire(fanInWireName)
            tempRepressor.addFanInWire(wire)
            fanOutWireToGateDict[fanInWireName].addFanOutWire(wire)
            allWires.append(wire)
        
        parIndex = findSmallestParentheses(circuitString)
        gateCount += 1
        wireCount += 1
    
    #Make the Output in a similar fashion to the repressors
    #Gate(name,Km,n,gateType,mB,pB,m_halfLife, p_halfLife)
    tempOut = Gate.Gate("Y",None,None,'Output',None,20,2,10)
    wire = Wire.Wire(circuitString)
    tempOut.addFanInWire(wire)
    fanOutWireToGateDict[circuitString].addFanOutWire(wire)
    Outputs.append(tempOut)
    allWires.append(wire)

    allGates = Inputs + Repressors + Outputs
    
    #Make the DAG with xMin of 0 and xMax of the period of the most significant
    #input. Arbitrarily chose 10 times xMax as the number of iterations to make
    #sure the graph will be fine enough.
    dag = DAG.DAG(0,initPeriod,10*initPeriod)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire) 
    
    return dag, allGates, Inputs, Repressors, Outputs
    
    
def makeGatesFromLibraries(Libraries,Inputs,Repressors,Outputs):
    """
    Allows you to choose which inputs, repressors and outputs to use from the
    set of libraries.
    """
    #Import the gate informations from the Libraries
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
    
    #Get a list of the valid choices for inputs, repressors, outputs
    inputsFromFileNames = []
    repressorsFromFileNames = []
    outputsFromFileNames = []
    for gateInfo in inputsFromFile:
        inputsFromFileNames.append(gateInfo["NAME"])
    for gateInfo in repressorsFromFile:
        repressorsFromFileNames.append(gateInfo["NAME"])
    for gateInfo in outputsFromFile:
        outputsFromFileNames.append(gateInfo["NAME"])
    
    #Get a string that can be used to display options
    inputOptions = ""
    repressorOptions = ""
    outputOptions = ""
    for i in range(len(inputsFromFileNames)):
        inputOptions = inputOptions + str(i) + " : " + inputsFromFileNames[i] + "\n"
    for i in range(len(repressorsFromFile)):
        repressorOptions = repressorOptions + str(i) + " : " + repressorsFromFileNames[i] + "\n"
    for i in range(len(outputsFromFileNames)):
        outputOptions = outputOptions + str(i) + " : " + outputsFromFileNames[i] + "\n"
    #Allow for the selection of inputs
    selectedInputs = []
    selectedRepressors = []
    selectedOutputs = []
    selectedInputNames = []
    selectedRepressorNames = []
    selectedOutputNames = []
    print "You must select",str(len(Inputs)),"unique inputs. Remember, order matters."
    while len(selectedInputs)<len(Inputs):
        print inputOptions
        index = raw_input("Type the number corresponding to the input you want: ")
        try:        
            index = int(index)
        except ValueError:
            print "This is not a valid selection"
            continue
        if index not in range(len(inputsFromFileNames)):
            print "This is not a valid selection"
        elif inputsFromFile[index] in selectedInputs:
            print "You have already chosen this option."
            print selectedInputNames
        else:
            selectedInputs.append(inputsFromFile[index])
            selectedInputNames.append(inputsFromFileNames[index])
            print ""


    print "You must select",str(len(Repressors)),"unique repressors. Remember, order matters."
    while len(selectedRepressors)<len(Repressors):
        print repressorOptions
        index = raw_input("Type the number corresponding to the repressor you want: ")
        try:        
            index = int(index)
        except ValueError:
            print "This is not a valid selection"
            continue
        if index not in range(len(repressorsFromFileNames)):
            print "This is not a valid selection"
        elif repressorsFromFile[index] in selectedRepressors:
            print "You have already chosen this option please select a different option that isn't in the following list."
            print selectedRepressorNames
        else:
            selectedRepressors.append(repressorsFromFile[index])
            selectedRepressorNames.append(repressorsFromFileNames[index])
            print ""
     
    
    print "You must select",str(len(Outputs)),"unique outputs. Remember, order matters."
    while len(selectedOutputs)<len(Outputs):
        print outputOptions
        index = raw_input("Type the number corresponding to the repressor you want: ")
        try:        
            index = int(index)
        except ValueError:
            print "This is not a valid selection"
            continue
        if index not in range(len(outputsFromFileNames)):
            print "This is not a valid selection"
        elif outputsFromFile[index] in selectedOutputs:
            print "You have already chosen this option."
            print selectedOutputNames
        else:
            selectedOutputs.append(outputsFromFile[index])
            selectedOutputNames.append(outputsFromFileNames[index])
            print ""
    
    #Make the gates
    print "Making gates from selections"
    for i in range(len(selectedInputs)):
        gateInfo = selectedInputs[i]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),gateInfo["Km"],gateInfo["n"],"Input",gateInfo["mB"],None,None,None)
        Inputsr.append(tempGate)
    for i in range(len(selectedRepressors)):
        gateInfo = selectedRepressors[i]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),gateInfo["Km"],gateInfo["n"],"Repressor",gateInfo["mB"],gateInfo["pB"],gateInfo["m_halfLife"],gateInfo["p_halfLife"])
        Repressorsr.append(tempGate)
    for i in range(len(selectedOutputs)):
        gateInfo = selectedOutputs[i]
        tempGate = Gate.Gate(str(gateInfo["NAME"]),None,None,"Output",None,gateInfo["pB"],gateInfo["m_halfLife"],gateInfo["p_halfLife"])
        Outputsr.append(tempGate)
    allGatesr = Inputsr + Repressorsr + Outputsr
    return allGatesr, Inputsr, Repressorsr, Outputsr
    
def performSwaps(dag,allGates,Inputs,allGatesr,Inputsr):
    """
    Given a dag, a set of gates in the dag and a new set of gates, it will swap
    in the new gates for the old gates. The number of each type of gate must be
    the same.
    """
    
    #Make swaps in dag
    print "Performing swaps"
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



def wrapper(circuitString, Libraries, fileLoc, makeBarGraph):
    """
    Takes in a string representation of a circuit, a set of libraries of inputs,
    repressors, and outputs, and a directory to store the intermediate json file.
    Prints the graphs generated for mRNA and protein concentrations generated
    by General.
    """
    #Make DAG from string
    dag, allGates, Inputs, Repressors, Outputs = makeDAGFromString(circuitString)
    #Get replacements for the gates
    allGatesr, Inputsr, Repressorsr, Outputsr = makeGatesFromLibraries(Libraries, Inputs, Repressors, Outputs)
    #Perform Swaps
    performSwaps(dag, allGates, Inputs, allGatesr, Inputsr)

    print dag
    #Write to Json File
    #fileLoc = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Example5.json'
    dag.writeToJson(fileLoc)
    #Make graphs from JsonFile
    General.generateDynamicCircuitGraphs(fileLoc, makeBarGraph)
    
repressilatorFileLoc = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan/repressilator.json"
placeToSaveRepressilator = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan/repressilatorDAG.json"
#
truthValueExampleFileLoc = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan\SplitByTruthValue_OR/01101001.json" #01101001
placeToSaveTruthValueExample = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\NetlistsFromBryan/01101001DAG.json"

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
         inputsDir = "C:\Users\Arinze\Documents\GitHub\MIT-BroadFoundry\dynamic_circuits\DynamicCircuits\Libraries/InputLibrary"+str(i)+".json"
         repressorsDir = "C:\Users\Arinze\Documents\GitHub\MIT-BroadFoundry\dynamic_circuits\DynamicCircuits\Libraries/RepressorLibrary"+str(i)+".json"
         outputsDir = "C:\Users\Arinze\Documents\GitHub\MIT-BroadFoundry\dynamic_circuits\DynamicCircuits\Libraries/OutputLibrary"+str(i)+".json"
         Libraries = [inputsDir, repressorsDir, outputsDir]
         wrapper(circuitString, Libraries, fileLoc, True)
         pause = raw_input("Finished "+str(i))
         
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
    
        