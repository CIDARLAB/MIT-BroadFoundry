# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 23:43:58 2015

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

defaultKm = 40
defaultn = 2
defaultmB = 30
defaultpB = 20

def makeNOTANORB():
    W1 = Wire.Wire('W1')
    W2 = Wire.Wire('W2')
    W3 = Wire.Wire('W3')
    W4 = Wire.Wire('W4')
    allWires = [W1,W2,W3,W4]
    a = Gate.Gate('a',40,2,'Input',30,None,None,None)
    b = Gate.Gate('b',40,2,'Input',30,None,None,None)
    G1 = Gate.Gate('G1',30000,2,'Repressor',30,20,2,10)
    G2 = Gate.Gate('G2',30000,2,'Repressor',30,20,2,10)
    y = Gate.Gate('y',None,None,'Output',None,20,2,10)
    allGates = [a,b,G1,G2,y]

    #Set Input value (fanOut)
    a.setFanOut([W1])
    a.setInputType(['inputs.squInput',1000,100,0,0])
    b.setFanOut([W2])
    b.setInputType(['inputs.squInput',500,100,0,0])
    #Set gate values (fanIn, fanOut, m_halfLife, p_halfLife)
    G1.addFanInWire(W1)
    G1.addFanInWire(W2)
    G1.setFanOut([W3])
    G2.addFanInWire(W3)
    G2.setFanOut([W4])
    #Set output values (fanIn, m_halfLife, p_halfLife)
    y.addFanInWire(W4)
    dag = DAG.DAG(0,2000,500)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
    print dag
    #print json.dumps(dag.formatJson(), sort_keys=True, indent=4, separators=(',', ': '))
    fileLoc = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Example3.json'
    dag.writeToJson(fileLoc)
    General.generateDynamicCircuitGraphs(fileLoc)
    
def practiceSwappingWires():
    #No Swap
    W1 = Wire.Wire('W1')
    W2 = Wire.Wire('W2')
    W3 = Wire.Wire('W3')
    allWires = [W1,W2,W3]

    G1 = Gate.Gate('G1',30000,2,'Repressor',30,20,2,10)
    G2 = Gate.Gate('G2',30000,2,'Repressor',30,20,2,10)
    G3 = Gate.Gate('G3',30000,2,'Repressor',30,20,2,10)
    allGates = [G1,G2,G3]
    #Set gate values (fanIn, fanOut, m_halfLife, p_halfLife)
    G1.addFanInWire(W1)
    G1.setFanOut([W2])
    G2.addFanInWire(W2)
    G2.setFanOut([W3])
    #Set output values (fanIn, m_halfLife, p_halfLife)
    dag = DAG.DAG(0,2000,500)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
    print G1.fanIn
    print G1.fanOut
    print G2.fanIn
    print G2.fanOut
    print G3.fanIn
    print G3.fanOut
    print dag
    
    
    #Remove a Wire
    W1 = Wire.Wire('W1')
    W2 = Wire.Wire('W2')
    W3 = Wire.Wire('W3')
    allWires = [W1,W2,W3]

    G1 = Gate.Gate('G1',30000,2,'Repressor',30,20,2,10)
    G2 = Gate.Gate('G2',30000,2,'Repressor',30,20,2,10)
    G3 = Gate.Gate('G3',30000,2,'Repressor',30,20,2,10)
    allGates = [G1,G2,G3]
    #Set gate values (fanIn, fanOut, m_halfLife, p_halfLife)
    G1.addFanInWire(W1)
    G1.setFanOut([W2])
    G1.setmRNAHalfLife(2)
    G1.setProteinHalfLife(10)
    G2.addFanInWire(W2)
    G2.setFanOut([W3])
    G2.setmRNAHalfLife(2)
    G2.setProteinHalfLife(10)
    #Remove the wires from G2
    G2.removeFanInWire(W2)
    G2.removeFanOutWire(W3)
    #Set output values (fanIn, m_halfLife, p_halfLife)
    dag = DAG.DAG(0,2000,500)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
    print G1.fanIn
    print G1.fanOut
    print G2.fanIn
    print G2.fanOut
    print G3.fanIn
    print G3.fanOut
    print dag
    
    
    #Swap a gate
    W1 = Wire.Wire('W1')
    W2 = Wire.Wire('W2')
    W3 = Wire.Wire('W3')
    allWires = [W1,W2,W3]

    G1 = Gate.Gate('G1',30000,2,'Repressor',30,20,2,10)
    G2 = Gate.Gate('G2',30000,2,'Repressor',30,20,2,10)
    G3 = Gate.Gate('G3',30000,2,'Repressor',30,20,2,10)
    allGates = [G1,G2,G3]
    #Set gate values (fanIn, fanOut, m_halfLife, p_halfLife)
    G1.addFanInWire(W1)
    G1.setFanOut([W2])
    G1.setmRNAHalfLife(2)
    G1.setProteinHalfLife(10)
    G2.addFanInWire(W2)
    G2.setFanOut([W3])
    G2.setmRNAHalfLife(2)
    G2.setProteinHalfLife(10)
    #Swap in G3 for G1
    G3.addFanInWire(W1)
    G3.addFanOutWire(W2)
    #Set output values (fanIn, m_halfLife, p_halfLife)
    dag = DAG.DAG(0,2000,500)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)
    print G1.fanIn
    print G1.fanOut
    print G2.fanIn
    print G2.fanOut
    print G3.fanIn
    print G3.fanOut
    print dag
    
def makeNOTANORBWithSwap():
    #Given a string representation of a circuit and the Input pattern.
    #Generate the DAG with 0 values for gate properties.
    W1 = Wire.Wire('W1')
    W2 = Wire.Wire('W2')
    W3 = Wire.Wire('W3')
    W4 = Wire.Wire('W4')
    allWires = [W1,W2,W3,W4]
    a = Gate.Gate('a',0,0,'Input',0,None,None,None)
    b = Gate.Gate('b',0,0,'Input',0,None,None,None)
    G1 = Gate.Gate('G1',0,0,'Repressor',0,0,0,0)
    G2 = Gate.Gate('G2',0,0,'Repressor',0,0,0,0)
    y = Gate.Gate('y',None,None,'Output',None,0,0,0)
    allGates = [a,b,G1,G2,y]
    Inputs = [a,b]
    #Set Input value (fanOut)
    a.setFanOut([W1])
    a.setInputType(['inputs.squInput',1000,100,0,0])    
    b.setFanOut([W2])
    b.setInputType(['inputs.squInput',500,100,0,0])
    #Set gate values (fanIn, fanOut, m_halfLife, p_halfLife)
    G1.addFanInWire(W1)
    G1.addFanInWire(W2)
    G1.setFanOut([W3])
    G2.addFanInWire(W3)
    G2.setFanOut([W4])
    #Set output values (fanIn, m_halfLife, p_halfLife)
    y.addFanInWire(W4)
    dag = DAG.DAG(0,2000,500)
    for gate in allGates:
        dag.addGate(gate)
    for wire in allWires:
        dag.addWire(wire)

    #Perform replaces
    #Make new gates from libraries
    ar = Gate.Gate('a',40,2,'Input',30,None,None,None)
    br = Gate.Gate('b',40,2,'Input',30,None,None,None)
    G1r = Gate.Gate('G1',30000,2,'Repressor',30,20,2,10)
    G2r = Gate.Gate('G2',30000,2,'Repressor',30,20,2,10)
    yr = Gate.Gate('y',None,None,'Output',None,20,2,10)
    allGatesr = [ar,br,G1r,G2r,yr]
    Inputsr = [ar,br]
    #Make swaps
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
    print dag
    #Write to Json File
    fileLoc = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ODE/Example4.json'
    dag.writeToJson(fileLoc)
    #Make graphs from JsonFile
    General.generateDynamicCircuitGraphs(fileLoc)