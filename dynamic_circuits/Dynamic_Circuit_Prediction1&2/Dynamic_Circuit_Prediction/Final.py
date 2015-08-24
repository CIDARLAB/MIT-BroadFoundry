# -*- coding: utf-8 -*-
import OptimalCircuit
import Optimize

repressilatorFileLoc = "JsonFiles/repressilator.json"
placeToSaveRepressilator = "JsonFiles/repressilatorGraph.json"
truthValueExampleFileLoc = "JsonFiles/SplitByTruthValue_OR/01101001.json" #01101001
exampleNetlist = ['NOR(W1,IN1,IN2)','NOR(W2,W1,IN3)','NOR(W3,IN1,IN3,IN4)','NOR(W4,W1,IN1,IN4)','NOR(W5,W2,W3,W4)']
placeToSaveExample = "JsonFiles/temp.json"
exampleCircuitString = "(((IN1.IN2).IN3).(IN1.IN3.IN4).((IN1.IN2).IN1.IN4))"
SR_Latch_FileLoc = "JsonFiles/SR_Latches/SRLatch.json"
SR_Latch_placeToSave = "JsonFiles/SR_Latches/SRLatch_Graph.json"
inputsDir = "JsonFiles/Libraries/InputLibrary2.json"
repressorsDir = "JsonFiles/Libraries/RepressorLibrary2.json"
outputsDir = "JsonFiles/Libraries/OutputLibrary2.json"
Libraries = [inputsDir, repressorsDir, outputsDir]



def makeGraphFromNetlist(netlist,placeToSave,makeBarGraph=True,makeGraphs=True):
    print OptimalCircuit.wrapperForNetlist(netlist,placeToSave,makeBarGraph=makeBarGraph,makeOtherGraphs=makeGraphs)

def makeGraphFromCircuitString(circuitString,placeToSave,makeBarGraph=True,makeGraphs=True):
    netlist = OptimalCircuit.convertCircuitStringToNetlist(circuitString)
    print OptimalCircuit.wrapperForNetlist(netlist,placeToSave,makeBarGraph=makeBarGraph,makeOtherGraphs=makeGraphs)

def optimizeNetlistWithLibraries(netlist, Libraries,smallestScoreAllowed=100,numTraj=20):
    answer = Optimize.findOptimalAssortmentHill(netlist,Libraries, smallestScoreAllowed,numTraj)
    print answer[0]
    print "Score =",answer[1]
    return answer
    
def optimizeCircuitStringWithLibraries(circuitString, Libraries,smallestScoreAllowed=100,numTraj=20):
    netlist = OptimalCircuit.convertCircuitStringToNetlist(circuitString)
    answer = Optimize.findOptimalAssortmentHill(netlist,Libraries, smallestScoreAllowed,numTraj)
    print answer[0]
    print "Score =",answer[1]
    return answer

def examplesForUse():
    makeGraphFromNetlist(repressilatorFileLoc,placeToSaveRepressilator)
    pause = raw_input("pause1")
    makeGraphFromNetlist(truthValueExampleFileLoc,placeToSaveExample)
    pause = raw_input("pause2")
    makeGraphFromNetlist(SR_Latch_FileLoc,SR_Latch_placeToSave)
    pause = raw_input("pause3")
    makeGraphFromNetlist(exampleNetlist,placeToSaveExample)
    pause = raw_input("pause4")
    makeGraphFromCircuitString(exampleCircuitString,placeToSaveExample)
    pause = raw_input("pause5")
    optimizeNetlistWithLibraries(truthValueExampleFileLoc, Libraries,smallestScoreAllowed=100,numTraj=5)
    
    
    