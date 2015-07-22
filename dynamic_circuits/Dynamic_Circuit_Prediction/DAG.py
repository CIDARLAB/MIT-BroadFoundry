# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 21:24:01 2015

@author: Arinze
"""

import numpy as np
import Gate
import Wire
import json
     
class DAG(object):
    
    def __init__(self,xMin,xMax,itr):
        self.inputs = []
        self.repressors = []
        self.outputs = []
        self.wires = []
        self.time_axis_params = {}
        self.time_axis_params['xMax'] = xMax
        self.time_axis_params['xMin'] = xMin
        self.time_axis_params['yMin'] = 'NONE'
        self.time_axis_params['yMax'] = 'NONE'
        self.time_axis_params['itr'] = itr
    
    def addWire(self,wire):
        assert type(wire) == Wire.Wire
        self.wires.append(wire)
    
    def removeWire(self,wire):
        assert type(wire) == Wire.Wire
        self.wires.remove(wire)
    
    def addGate(self,gate):
        assert type(gate)==Gate.Gate
        if gate.getGateType() == 'Input':
            self.inputs.append(gate)
        elif gate.getGateType() == 'Repressor':
            self.repressors.append(gate)
        elif gate.getGateType() == 'Output':
            self.outputs.append(gate)
        
    
    def removeGate(self,gate):
        assert type(gate)==Gate.Gate
        if gate.getGateType() == 'Input':
            self.inputs.remove(gate)
        elif gate.getGateType() == 'Repressor':
            self.repressors.remove(gate)
        elif gate.getGateType() == 'Output':
            self.outputs.remove(gate)
            
    def getAllGates(self):
        gates = self.inputs + self.repressors + self.outputs
        return gates
    
    def getAllWires(self):
        return self.wires
        
    def formatJson(self):
        listVals = [self.time_axis_params]
        for gate in self.inputs:
            listVals.append(gate.formatJson())
        for gate in self.repressors:
            listVals.append(gate.formatJson())
        for gate in self.outputs:
            listVals.append(gate.formatJson())
        return listVals
    def writeToJson(self,fileLoc):
        listVals = self.formatJson()
        myFile = open(fileLoc,'w')
        json.dump(listVals, myFile, sort_keys=True, indent=4, separators=(',', ': '))
        myFile.close()
    def gateCompare(self,gate1,gate2):
        """
        Used to sort gates assuming the name of the gate is a single letter
        followed by the gate number.
        """
                
        if gate1.getDist()>gate2.getDist():
            return 1
        elif gate1.getDist()==gate2.getDist():
            return 0
        else: #gate1Num<gate2Num
            return -1

    def sortGates(self):
        self.repressors
        sortedRepressors = []
        while len(sortedRepressors) != len(self.repressors):
            for rep in self.repressors:
                if len(rep.getFanIn())==2:
                    if rep not in sortedRepressors and \
                    (rep.getFanIn()[0].getFrom() in sortedRepressors \
                    or rep.getFanIn()[0].getFrom().getGateType()=="Input") and\
                    (rep.getFanIn()[1].getFrom() in sortedRepressors or \
                    rep.getFanIn()[1].getFrom().getGateType()=="Input"):
                        sortedRepressors.append(rep)
                if len(rep.getFanIn())==1:
                    if rep not in sortedRepressors and \
                    (rep.getFanIn()[0].getFrom() in sortedRepressors or \
                    rep.getFanIn()[0].getFrom().getGateType()=="Input"):
                        sortedRepressors.append(rep)
        self.repressors = sortedRepressors
    def setGateDist(self):
        """
        Assumes you have a fully connected graph and assigns distances from the
        inputs
        """
        #Visited property is to prevent infinite loops with sequential
        #circuits
        for gate in self.repressors:
            gate.setVisited(False)
        for gate in self.outputs:
            gate.setVisited(False)
        for gate in self.outputs:
            self.recursivelyFindInputs(gate)
        
    def recursivelyFindInputs(self,gate):
        if gate.wasVisited():
            return
        elif gate.getGateType()=="Input":
            return
        else:
            gate.setVisited(True)
            vals = []
            for wire in gate.getFanIn():
                self.recursivelyFindInputs(wire.getFrom())
                vals.append(wire.getFrom().getDist())
            gate.setDist(max(vals)+1)
            
        
    def __str__(self):
        result = 'Inputs:\n'
        for gate in self.inputs:
            result = result + gate.getName() + '\n'
        result = result + 'Repressors:\n'
        for gate in self.repressors:
            result = result + gate.getName() + '\n'
        result = result + 'Outputs:\n'
        for gate in self.outputs:
            result = result + gate.getName() + '\n\n'
            
        result = result + 'NETLIST:\n'
        self.setGateDist()
        self.repressors.sort(self.gateCompare)
        self.outputs.sort(self.gateCompare)
#        self.sortGates()
        for gate in self.repressors:
            gateStr = ""
            if len(gate.getFanIn()) == 2:
                gateStr = gateStr + "NOR("
            if len(gate.getFanIn()) == 1:
                gateStr = gateStr + "NOT("
            gateStr = gateStr + gate.getName() + ","
            for w in gate.getFanIn():
                gateStr = gateStr + w.getFrom().getName() + ","
            gateStr = gateStr[:-1] + "); distance = "+str(gate.getDist())+"\n"
            result = result + gateStr
        for gate in self.outputs:
            gateStr = ""
            if len(gate.getFanIn()) == 2:
                gateStr = gateStr + "OR("
            if len(gate.getFanIn()) == 1:
                gateStr = gateStr + "BUF("
            gateStr = gateStr + gate.getName() + ","
            for w in gate.getFanIn():
                gateStr = gateStr + w.getFrom().getName() + ","
            gateStr = gateStr[:-1] + "); distance = "+str(gate.getDist())+"\n"
            result = result + gateStr
        return result
