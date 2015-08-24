# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 21:24:01 2015

@author: Arinze
"""

import numpy as np
import Gate
import Wire
import json
import copy
     
class DAG(object):
    
    def __init__(self,xMin,xMax,itr):
        self.inputs = []
        self.intermediates = []
        self.outputs = []
        self.wires = []
        self.time_axis_params = {}
        self.time_axis_params['xMax'] = xMax
        self.time_axis_params['xMin'] = xMin
        self.time_axis_params['yMin'] = 'NONE'
        self.time_axis_params['yMax'] = 'NONE'
        self.time_axis_params['itr'] = itr
        self.isSequential = False
    
    def addWire(self,wire):
        """
        adds a wire to the DAG if it is not already in the DAG.
        """
        assert type(wire) == Wire.Wire
        if wire not in self.wires:
            self.wires.append(wire)
        else:
            print "This wire is already in the DAG"
    
    def removeWire(self,wire):
        """
        removes a wire from the DAG if it is in the DAG
        """
        assert type(wire) == Wire.Wire
        if wire in self.wires:
            self.wires.remove(wire)
        else:
            print "This wire is not in the DAG"
    
    def addGate(self,gate):
        """
        adds a gate to the DAG if it is not already in the DAG
        """
        #This format is used in case the name of the gate is changed after 
        #addition to the DAG
        assert type(gate)==Gate.Gate
        allGates = self.inputs+self.intermediates+self.outputs
        if gate not in allGates:
            if gate.getGateType() == 'Input':
                self.inputs.append(gate)
            elif gate.getGateType() == 'Repressor' or gate.getGateType() == 'Activator':
                self.intermediates.append(gate)
            elif gate.getGateType() == 'Output':
                self.outputs.append(gate)
        else:
            if gate in self.inputs:
                print "This gate is already present in the DAG inputs."
            elif gate in self.intermediates:
                print "This gate is already present in the DAG intermediates."
            elif gate in self.outputs:
                print "This gate is already present in the DAG outputs."
    
    def removeGate(self,gate):
        """
        removes a gate from the DAG if it is in the DAG
        """
        #This format is used in case the name of the gate is changed after 
        #addition to the DAG
        assert type(gate)==Gate.Gate
        try:
            self.inputs.remove(gate)
        except ValueError:
            try:
                self.intermediates.remove(gate)
            except ValueError:
                try:
                    self.outputs.remove(gate)
                except ValueError:
                    print "This gate is not in the DAG"
            
    def getAllGates(self):
        """
        Returns a list of all gates in the circuit
        """
        gates = self.inputs + self.intermediates + self.outputs
        return gates
    
    def getAllWires(self):
        """
        Returns a list of all wires in the DAG
        """
        return self.wires
    def hasLoop(self):
        """
        Returns true if the DAG is sequential
        """
        return self.isSequential
        
    def swapGates(self,gate1,gate2):
        """
        Given two gates, performs a swap. Setting the fanIn for one as the 
        fanIn for the other and vice versa and the same for fanOut.
        """
        assert type(gate1)==Gate.Gate and type(gate2)==Gate.Gate
        fanIn1 = copy.copy(gate1.getFanIn())
        fanOut1 = copy.copy(gate1.getFanOut())
        fanIn2 = copy.copy(gate2.getFanIn())
        fanOut2 = copy.copy(gate2.getFanOut())
        
        gate1.setFanIn(fanIn2)
        gate1.setFanOut(fanOut2)
        gate2.setFanIn(fanIn1)
        gate2.setFanOut(fanOut1)
        
        allGates = self.inputs + self.intermediates + self.outputs
        gate1InDAG = False
        gate2InDAG = False
        if gate1 in allGates:
            gate1InDAG = True
        if gate2 in allGates:
            gate2InDAG = True
        
        if gate1InDAG == gate2InDAG:
            return
        elif not gate1InDAG and gate2InDAG:
            self.addGate(gate1)
            self.removeGate(gate2)            
            
        elif gate1InDAG and not gate2InDAG:
            self.addGate(gate2)
            self.removeGate(gate1)
        
    def formatJson(self):
        """
        formats all of the gates so as to be easily exported to a json file
        """
        self.setGateDist()
        self.intermediates.sort(self.gateCompare)
        self.outputs.sort(self.gateCompare)
        self.setGateTruthValues()
        listVals = [self.time_axis_params]
        for gate in self.inputs:
            listVals.append(gate.formatJson())
        for gate in self.intermediates:
            listVals.append(gate.formatJson())
        for gate in self.outputs:
            listVals.append(gate.formatJson())
        return listVals
        
    def writeToJson(self,fileLoc):
        """
        writes the contents of the dag to a json with pretty print in a format
        that can be recognized by General.py
        """
        listVals = self.formatJson()
        myFile = open(fileLoc,'w')
        json.dump(listVals, myFile, sort_keys=True, indent=4, separators=(',', ': '))
        myFile.close()
        
    def gateCompare(self,gate1,gate2):
        """
        Used to sort gates based on their distance from the furthest input 
        discounting loops
        """
                
        if gate1.getDist()>gate2.getDist():
            return 1
        elif gate1.getDist()==gate2.getDist():
            return 0
        else: #gate1Num<gate2Num
            return -1

    def setGateDist(self):
        """
        Assumes you have a fully connected graph and assigns distances from the
        inputs
        """
        #Visited property is to prevent infinite loops with sequential
        #circuits
        for gate in self.intermediates:
            gate.setVisited(False)
        for gate in self.outputs:
            gate.setVisited(False)
        for gate in self.intermediates:
            self.recursivelyFindInputs(gate)
        for gate in self.outputs:
            self.recursivelyFindInputs(gate)
        
    def recursivelyFindInputs(self,gate):
        """
        a helper function that does the actual finding of the distance to the 
        furthest input discounting loops
        """
        #Stops loops
        if gate.wasVisited():
            return
        #Stops at inputs
        elif gate.getGateType()=="Input":
            return
        #save the distance from each of the gates going into this
        #take the max of those values and add 1
        else:
            gate.setVisited(True)
            vals = []
            for wire in gate.getFanIn():
                self.recursivelyFindInputs(wire.getFrom())
                vals.append(wire.getFrom().getDist())
            gate.setDist(max(vals)+1)
            
    def setGateTruthValues(self):
        """
        Assumes you have a fully connected graph and assigns truthValues to 
        the gates in a recursive way similar to finding the distance to the 
        nearest input
        """
        #Visited property is to prevent infinite loops with sequential
        #circuits
        for gate in self.intermediates:
            gate.setVisited(False)
        for gate in self.outputs:
            gate.setVisited(False)
            
        for gate in self.intermediates:
            gate.setExpectedProtein(None)
        for gate in self.outputs:
            gate.setExpectedProtein(None)
        for gate in self.intermediates:
            self.recursivelyFindInputsAndSetTruthValues(gate)
        for gate in self.outputs:
            self.recursivelyFindInputsAndSetTruthValues(gate)
        
    def recursivelyFindInputsAndSetTruthValues(self,gate):
        """
        This function is what actually does the truth value assignment
        """

        if gate.wasVisited():
            if gate.getExpectedProtein()==None:
                self.isSequential = True
            return
        elif gate.getExpectedProtein()!=None and gate.getExpectedPromoter()!=None:
            gate.setVisited(True)
            return
        else:
            gate.setVisited(True)
            vals = []
            for wire in gate.getFanIn():
                self.recursivelyFindInputsAndSetTruthValues(wire.getFrom())
                vals.append(wire.getFrom().getExpectedPromoter())
            #If there is a loop, set all truth values dependent on that loop 
            #to None otheriwse find the truthValue for the protein
            try:
                truthValue = []
                numFanIn = len(vals)
                for i in range(len(vals[0])):
                    isOn = False
                    for j in range(numFanIn):
                        if vals[j][i]=="1":
                            isOn = True
                            break
                    if isOn:
                        truthValue.append("1")
                    elif not isOn:
                        truthValue.append("0")
                gate.setExpectedProtein(truthValue)
            except TypeError:
                gate.setExpectedProtein(None)
                self.isSequential = True

    def __str__(self):
        result = 'Inputs:\n'
        for gate in self.inputs:
            result = result + gate.getName() + '\n'
        result = result + 'Intermediates:\n'
        for gate in self.intermediates:
            result = result + gate.getName() + '\n'
        result = result + 'Outputs:\n'
        for gate in self.outputs:
            result = result + gate.getName() + '\n\n'
            
        result = result + 'NETLIST:\n'
        self.setGateDist()
        self.intermediates.sort(self.gateCompare)
        self.outputs.sort(self.gateCompare)
        self.setGateTruthValues()
        for gate in self.intermediates:
            result = result + str(gate) + "\n"
        for gate in self.outputs:
            result = result + str(gate) + "\n"
        return result
