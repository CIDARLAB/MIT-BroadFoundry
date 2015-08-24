# -*- coding: utf-8 -*-
import LoadNetlist
import Objects
import LogicOperations

#Latches
Gated_D_Latch_FileLoc = "JsonFiles/SR_Latches/Gated_D_Latch.json"
Gated_D_Latch_noCInv_FileLoc = "JsonFiles/SR_Latches/Gated_D_Latch_noCInv.json"
Negative_Edge_Triggered_D_Flip_Flop_FileLoc = "JsonFiles/SR_Latches/Negative_Edge_Triggered_D_Flip_Flop.json"
Transparent_D_Latch_noDinv_FileLoc = "JsonFiles/SR_Latches/Transparent_D_Latch_noDinv.json"
Transparent_D_Latch_noInpinv_FileLoc = "JsonFiles/SR_Latches/Transparent_D_Latch_noInpinv.json"
#TestCase
Test1_FileLoc = "JsonFiles/SR_Latches/Test1.json"
Test2_FileLoc = "JsonFiles/SR_Latches/Test2.json"

exampleNetlist = ['NOR(w1,a,b)','NOR(w2,w1,c)','NOT(w3,w2)','OR(w4,w3,w2)','BUF(y,w4)']

def checkStates(netlistDir):
    graph = LoadNetlist.makeGraphFromNetlist(netlistDir)[0]
#    graph.setGateStates([1,0,0,0])
    graph.runUntilStableOrRepeat(0)
    
    print graph
    