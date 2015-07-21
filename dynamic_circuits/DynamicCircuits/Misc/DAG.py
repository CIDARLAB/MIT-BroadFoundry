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
        return self.gates
    
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
    def __str__(self):
        result = 'Inputs:\n'
        for gate in self.inputs:
            result = result + gate.getName() + '\n'
        result = result + 'Repressors:\n'
        for gate in self.repressors:
            result = result + gate.getName() + '\n'
        result = result + 'Outputs:\n'
        for gate in self.outputs:
            result = result + gate.getName() + '\n'
        result = result + 'Wires:\n'
        for wire in self.wires:
            result = result + str(wire) + '\n'
        return result
