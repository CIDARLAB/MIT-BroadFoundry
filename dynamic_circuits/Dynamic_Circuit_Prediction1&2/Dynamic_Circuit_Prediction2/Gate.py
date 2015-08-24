# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 22:30:10 2015

@author: Arinze
Edited by Alex on 7/30/2015
"""
import numpy as np
import Wire
import Graph
import inputs
import copy

class Gate(object):
    def __init__(self,name,gateType,Km,n,Pmin,Pmax,halfLife,REUconv=None):#,mB,pB,m_halfLife,p_halfLife):
        assert type(name) == str
        assert gateType in ['Repressor','Input','Output','Activator']
        #Name the gate
        self.name = name
        #Default the fanIn and fanOut to None until they are set
        self.fanIn = []
        self.fanOut = []
        #Repressor, Input, Output
        self.gateType = gateType
        #binding efficiency of the protein produced to the promoter of subsequent gene
        self.Km = None
        #hill coefficient for the protein binding to the promoter of subsequent gene
        self.n = None
        self.Pmin = None
        self.Pmax = None
        try:
            self.REUconv = float(REUconv)
        except TypeError:
            self.REUconv = REUconv
        try:
            self.halfLife = float(halfLife)
        except TypeError:
            self.halfLife = halfLife
            
        self.degRate = None
        self.REUi = None
        #To adjust the input equation for Inputs
        self.inputType = None
        #Properties to be manipulated by Graph
        self.dist = 0    
        self.visited = True
        self.expectedProtein = None
        self.expectedPromoter = None
        
        if self.gateType != 'Output':
            self.Km = float(Km)
            self.n = float(n)
            self.Pmin = float(Pmin)
            self.Pmax = float(Pmax)
#            self.mB = None
        if self.gateType != 'Input':
            self.degRate = np.log(2)/self.halfLife
            self.REUi = 0.0
            self.visited = False
                
    #Setter functions
    def setName(self,name):
        '''
        Takes in a string that will be set to the name of the wire.
        Allows you to change the name
        '''
        #Make sure input type is correct
        assert type(name) == str
        #Set value
        self.name = name
        
    def setFanIn(self,fanIn):
        '''
        Allows you to reset the fanIn to either an empty list or a list of wires
        '''
        #Inputs cannot have a fanIn
        if self.gateType == "Input":
            print "An input gate cannot have a fanIn"
            return
        #Make sure the input type is correct
        assert type(fanIn) == list
        for wire in fanIn:
            assert type(wire) == Wire.Wire
        fanIn = copy.copy(fanIn)
        #Reset the To for the wires that were already in the fanIn
        for wire in self.fanIn:
            wire.setTo(None)
        #Set the values
        self.fanIn = []
        while len(fanIn)!=0:
            # Takes the 1st wire from list fanIn and removes its connection 
            # to this gate. Upon removal from the list, the next wire in the 
            # list is placed into index 0 of the list.
            wire = fanIn[0]
            #Let the wire know it is connected to the gate
            wire.setTo(self)
            try:
                fanIn.remove(wire)
            except ValueError:
                pass
            #Let the gate know it is connected to the wire
            self.fanIn.append(wire)
        
    def addFanInWire(self,wire):
        '''
        Adds a Wire to the fanIn
        '''
        if self.gateType == "Input":
            print "An Input cannot have a fanIn"
            return
        #Make sure the input type is correct
        assert type(wire) == Wire.Wire
        #Set the values
        wire.setTo(self)
        self.fanIn.append(wire)
    
    def removeFanInWire(self,wire):
        '''
        Removes a Wire from fanIn. Not to be used by itself. Only to be called
        by wires.
        '''
        assert type(wire) == Wire.Wire
        self.fanIn.remove(wire)   
        wire.setTo(None)
        
    def addFanOutWire(self,wire):
        '''
        Add a Wire to the fanOut
        '''
        if self.gateType == "Output":
            print "An Output gate cannot have a fanOut"
            return
        #Make sure the input type is correct
        assert type(wire) == Wire.Wire
        #Set the values
        wire.setFrom(self)
        self.fanOut.append(wire)
        
    def removeFanOutWire(self,wire):
        '''
        Removes a Wire from fanOut. Not to be used by itself, only to be called
        by wires.
        '''
        assert type(wire) == Wire.Wire
        self.fanOut.remove(wire)
        wire.setFrom(None)
        
    def setFanOut(self,fanOut):
        '''
        Allows you to reset the fanOut to either an empty list or a list of wires
        '''
        if self.gateType == "Output":
            print "An Output cannot have a fanOut"
            return
        #Make sure the input type is correct
        assert type(fanOut) == list
        for wire in fanOut:
            assert type(wire) == Wire.Wire
        fanOut = copy.copy(fanOut)
        #Reset the To for the wires that were already in the fanOut
        for wire in self.fanOut:
            wire.setFrom(None)
        #set the values
        self.fanOut = []
        #This format is used in case fanOut is being mutated as we iterate 
        #through it. Instead of risking possible mutation, we purposefully
        #mutate it if it isnt being mutated.
        while len(fanOut)!=0:
            wire = fanOut[0]
            wire.setFrom(self)
            try:
                fanOut.remove(wire)
            except ValueError:
                pass
            self.fanOut.append(wire)
        
    def setKm(self,Km):
        '''
        To Change the Km value
        '''
        self.Km = float(Km)
        
    def setn(self,n):
        '''
        To Change the n value
        '''
        self.n = float(n)
    def setPmin(self,Pmin):
        '''
        To Change the Pmin value
        '''
        self.Pmin = float(Pmin)
        
    def setPmax(self,Pmax):
        '''
        To Change the Pmax value
        '''
        self.Pmax = float(Pmax)
        
    def setHalfLife(self,halfLife):
        '''
        Sets the half-life and updates the degradation rate
        '''
        try:
            self.halfLife = float(halfLife)
            self.degRate = np.log(2)/self.halfLife
        except TypeError:
            print "This is not a valid half life"
            return

    def setREUConv(self,REUconv):
        '''
        To Change the Pmax value
        '''
        self.REUconv = float(REUconv)
    
    def setInputType(self, inputType):
        '''
        Sets the input to an equation from inputs
        '''
        if inputType==None:
            self.inputType = inputType
            return
        assert self.gateType=='Input' and type(inputType)==list
        self.availableInputTypes = ['inputs.sinInput', 'inputs.linInput','inputs.sawInput','inputs.squInput', 'inputs.stepFunction', 'inputs.stepInput','inputs.specInput',\
        'sinInput', 'linInput','sawInput','squInput', 'stepFunction', 'stepInput','specInput']
        self.numInputs = [4,2,3,4,3,4,4,4,2,3,4,3,4,4]
        assert inputType[0] in self.availableInputTypes
        #Make sure we have the correct number of inputs
        i = self.availableInputTypes.index(inputType[0])
        assert len(inputType) <= self.numInputs[i]+1
        
        for i in range(1,len(inputType)):
            if inputType[0] in ['specInput','inputs.specInput'] and i==len(inputType)-1:
                continue
            inputType[i] = float(inputType[i])
        self.inputType = inputType
        
    def setREUi(self, REUi):
        assert self.gateType!='Input'
        self.REUi = REUi
        
    def setDist(self, d):
        self.dist = int(d)
        
    def setVisited(self, visited):
        assert type(visited)==bool
        self.visited = visited

    def setExpectedProtein(self,tv):
        """
        Sets the expected protein and promoter truth values. Promoter truth 
        values are the opposite of the protein truth value if the gene is a
        repressor. Outputs do not have a promoter.
        """
        assert type(tv) == list or tv == None
        self.expectedProtein = tv
        if tv!=None:
            if self.gateType == "Input" or self.gateType == "Activator":
                self.setExpectedPromoter(tv)
            elif self.gateType == "Repressor":
                self.setExpectedPromoter(self.invertTruth(tv))
            elif self.gateType == "Output":
                self.setExpectedPromoter(None)
        elif tv==None:
            self.setExpectedPromoter(tv)
            
        
    def setExpectedPromoter(self,tv):
        """
        Set the expected setting of the promoter that comes after the gene as 
        a truth value.
        """
        assert type(tv) == list or tv == None
        self.expectedPromoter = tv
        
    def invertTruth(self, tv):
        """
        Used to invert truth values. Takes in a list of '1' and '0's and 
        inverts the value at each position. Does not alter any of the 
        properties of the gate.
        """
        newTruthValue = []
        for i in tv:
            if i == "0":
                newTruthValue.append("1")
            elif i == "1":
                newTruthValue.append("0")
        return newTruthValue

    def getName(self):
        return self.name
    def getFanIn(self):
        return self.fanIn
    def getFanOut(self):
        return self.fanOut
    def getKm(self):
        return self.Km
    def getn(self):
        return self.n
    def getGateType(self):
        return self.gateType
    def getPmax(self):
        return self.Pmax
    def getPmin(self):
        return self.Pmin
    def getHalfLife(self):
        return self.halfLife
    def getDegRate(self):
        return self.degRate
    def getREUconv(self):
        return self.REUconv
    def getREUi(self):
        return self.REUi        
    def getInputType(self):
        return self.inputType
    def getDist(self):
        return self.dist
    def wasVisited(self):
        return self.visited
    def getExpectedProtein(self):
        return self.expectedProtein
    def getExpectedPromoter(self):
        return self.expectedPromoter

    def formatJson(self):
        """
        creates a dictionary of all important properties necessary for
        General.py to produce graphs for easier export to a json file
        """
        gateDescription = {}
        gateDescription['NAME'] = self.name
        gateDescription['DIST'] = self.dist
        gateDescription['EXPECTED_PROTEIN'] = self.expectedProtein
        gateDescription['EXPECTED_PROMOTER'] = self.expectedPromoter
        gateDescription['REUconv'] = self.REUconv
        gateDescription['REUi'] = self.REUi
        gateDescription['degRate'] = self.degRate
        gateDescription['halfLife'] = self.halfLife
        gateDescription['Pmin'] = self.Pmin
        gateDescription['Pmax'] = self.Pmax
        gateDescription['Km'] = self.Km
        gateDescription['n'] = self.n
        if self.gateType=='Input':
            gateDescription['TYPE'] = 'INPUT'
            gateDescription['INPUT'] = self.inputType
            gateDescription['INPUT_EFFECT'] = "ACTIVATE"
        else:
            gateDescription['INPUT'] = []
            for i in range(len(self.fanIn)):
                gateDescription['INPUT'].append(self.fanIn[i].getFrom().getName())
            
            if self.gateType=="Input" or self.gateType=="Activator":
                gateDescription['INPUT_EFFECT'] = "ACTIVATE"
            elif self.gateType=="Repressor":
                gateDescription['INPUT_EFFECT'] = "REPRESS"
            if len(self.fanIn)==1:
                if self.gateType=='Repressor':
                    gateDescription['TYPE'] = 'NOT'
                elif self.gateType=='Activator' or self.gateType=='Output':
                    gateDescription['TYPE'] = 'BUFFER'
            elif len(self.fanIn)>=2:
                if self.gateType=='Repressor':
                    gateDescription['TYPE'] = 'NOR'
                elif self.gateType=='Activator' or self.gateType=='Output':
                    gateDescription['TYPE'] = 'OR'
        return gateDescription
        
    def __repr__(self):
        return self.name
    
    def __str__(self):
        """
        prints the gate in the way it would appear in a netlist
        """
        gateStr = ""
        if self.gateType == "Input":
            gateStr = gateStr + "INPUT("
        elif len(self.fanIn) >= 2 and self.gateType == "Repressor":
            gateStr = gateStr + "NOR("
        elif len(self.fanIn) == 1 and self.gateType == "Repressor":
            gateStr = gateStr + "NOT("
        elif len(self.fanIn) >= 2 and (self.gateType == "Output" or self.gateType == "Activator"):
            gateStr = gateStr + "OR("
        elif len(self.fanIn) == 1 and (self.gateType == "Output" or self.gateType == "Activator"):
            gateStr = gateStr + "BUF("
        else:
            gateStr = gateStr + "???("
        gateStr = gateStr + self.name + ","
        for w in self.fanIn:
            try:
                gateStr = gateStr + w.getFrom().getName() + ","
            except TypeError:
                gateStr = gateStr + "None" + ","
        #Remove the ending comma"
        gateStr = gateStr[:-1] + ")"
        return gateStr
  