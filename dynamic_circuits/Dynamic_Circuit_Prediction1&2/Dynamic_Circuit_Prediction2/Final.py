# -*- coding: utf-8 -*-
import OptimalCircuit
import Optimize

repressilatorFileLoc = "JsonFiles/repressilator.json"

truthValueExampleFileLoc = "JsonFiles/SplitByTruthValue_OR/01101001.json" #01101001
truthValueExampleFileLoc2 = "JsonFiles/SplitByTruthValue_OR/11110010.json" #11110010

#exampleNetlist = ['NOR(W1,IN1,IN2)','NOR(W2,W1,IN3)','NOR(W3,IN1,IN3,IN4)','NOR(W4,W1,IN1,IN4)','NOR(W5,W2,W3,W4)']
exampleNetlist = ['NOT(W1,IN1)', 'NOT(W2,IN2)', 'NOR(W3,W2,W1)']

exampleInputsList = ["Ara","IPTG"]
exampleIntermediatesList = ["E1_BetI","B1_BM3R1","A1_AmtR"]
exampleOutputList = ["BFP"]
genesList = [exampleInputsList,exampleIntermediatesList,exampleOutputList]

exampleCircuitString = "(((IN1.IN2).IN3).(IN1.IN3.IN4).((IN1.IN2).IN1.IN4))"

SR_Latch_FileLoc = "JsonFiles/SR_Latches/SRLatch.json"

inputsDir = "JsonFiles/Libraries/InputLibrary.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary.json"
Libraries = [inputsDir, repressorsDir, outputsDir]

#placeToSave converts a netlist (list of strings) to a JSON file.
#'General.py' reads the json file and makes a Graph.
def makeGraphFromNetlist(netlist,makeBarGraph=True,makeGraphs=True,useDefaultWaveForm=True,genesToUse=None):
    placeToSave = "JsonFiles/singleGraphTempFile1.json"
    print OptimalCircuit.wrapperForNetlist(netlist,placeToSave,makeBarGraph=makeBarGraph,makeOtherGraphs=makeGraphs,useDefaultInput=useDefaultWaveForm,genesToUse=genesToUse)

def makeGraphFromCircuitString(circuitString,makeBarGraph=True,makeGraphs=True,useDefaultWaveForm=True,genesToUse=None):
    netlist = OptimalCircuit.convertCircuitStringToNetlist(circuitString)
    placeToSave = "JsonFiles/singleGraphTempFile2.json"
    print OptimalCircuit.wrapperForNetlist(netlist,placeToSave,makeBarGraph=makeBarGraph,makeOtherGraphs=makeGraphs,useDefaultInput=useDefaultWaveForm,genesToUse=genesToUse)

#combinational only, not sequential
def optimizeNetlistWithLibraries(netlist,Libraries=Libraries,smallestScoreAllowed=10,numTraj=5):
    answer = Optimize.findOptimalAssortmentHill(netlist,Libraries, smallestScoreAllowed,numTraj)
    print answer[0]
    print "Score =",answer[1]
    return answer
    
#combinational only, not sequential
def optimizeCircuitStringWithLibraries(circuitString,Libraries=Libraries,smallestScoreAllowed=10,numTraj=5):
    netlist = OptimalCircuit.convertCircuitStringToNetlist(circuitString)
    answer = Optimize.findOptimalAssortmentHill(netlist,Libraries, smallestScoreAllowed,numTraj)
    print answer[0]
    print "Score =",answer[1]
    return answer

#combinational only, not sequential
def optimizeNetlistWithLibrariesTimed(netlist,Libraries=Libraries,smallestScoreAllowed=10,numTraj=5,maxTime=300):
    answer = Optimize.findOptimalAssortmentHillTimed(netlist,Libraries, smallestScoreAllowed,numTraj,maxTime)
    print answer[0]
    print "Score =",answer[1]
    return answer
    
#combinational only, not sequential
def optimizeCircuitStringWithLibrariesTimed(circuitString,Libraries=Libraries,smallestScoreAllowed=10,numTraj=5,maxTime=300):
    netlist = OptimalCircuit.convertCircuitStringToNetlist(circuitString)
    answer = Optimize.findOptimalAssortmentHillTimed(netlist,Libraries, smallestScoreAllowed,numTraj,maxTime)
    print answer[0]
    print "Score =",answer[1]
    return answer    
    
def examplesForUse():
    makeGraphFromNetlist(repressilatorFileLoc)
    pause = raw_input("pause1")
    makeGraphFromNetlist(truthValueExampleFileLoc)
    pause = raw_input("pause2")
    makeGraphFromNetlist(SR_Latch_FileLoc)
    pause = raw_input("pause3")
    makeGraphFromCircuitString(exampleCircuitString)
    pause = raw_input("pause4")
    makeGraphFromNetlist(exampleNetlist)
    pause = raw_input("pause5")
    makeGraphFromNetlist(exampleNetlist,genesToUse=genesToUse)
    pause = raw_input("pause6")
    optimizeNetlistWithLibraries(truthValueExampleFileLoc2,smallestScoreAllowed=3)