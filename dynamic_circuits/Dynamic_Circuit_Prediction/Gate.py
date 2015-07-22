# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 22:30:10 2015

@author: Arinze
"""
import numpy as np
import Wire
import DAG
import inputs



class Gate(object):
    def __init__(self,name,Km,n,gateType,mB,pB,m_halfLife, p_halfLife):
        assert type(name) == str
        assert gateType == 'Repressor' or gateType == 'Input' or gateType == 'Output'
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
        #mRNA induction coefficient for the promoter of subsequent gene
        self.mB = None
        #protein induction coefficient for this gene
        self.pB = None
        #half life of the mRNA produced by self
        self.m_halfLife = None
        #half life of the protein produced by self
        self.p_halfLife = None
        #The degradation rate of the mRNA produced by this gene
        self.am = None
        #The degradation rate of the protein produced by this gene
        self.ap = None
        #initial mRNA and protein concentrations
        self.Mi = None
        self.Pi = None
        #To adjust the input equation for Inputs
        self.inputType = None
        self.dist = 0    
        self.visited = True
        self.expectedProtein = None
        self.expectedPromoter = None
        
        if self.gateType != 'Output':
            self.Km = float(Km)
            self.n = float(n)
            self.mB = float(mB)
        if self.gateType != 'Input':
            self.pB = float(pB)
            self.m_halfLife = float(m_halfLife)
            self.p_halfLife = float(p_halfLife)
            self.am = np.log(2)/self.m_halfLife
            self.ap = np.log(2)/self.p_halfLife
            self.Mi = 0.0
            self.Pi = 0.0
            #self.dist = 9999 #arbitrary high value
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
        assert type(fanIn) == list and len(fanIn)<=2
        for wire in fanIn:
            assert type(wire) == Wire.Wire
        #Reset the To for the wires that were already in the fanIn
        for wire in self.fanIn:
            wire.setTo(None)
        #set the values
        self.fanIn = []
        while len(fanIn)!=0:
            wire = fanIn[0]
            wire.setTo(self)
            self.fanIn.append(wire)
        
    def addFanInWire(self,wire):
        '''
        Add a Wire to the fanIn
        '''
        if self.gateType == "Input":
            print "An input gate cannot have a fanIn"
            return
        if len(self.fanIn)==2:
            print "The fanIn already has the maximum number of wires."
            return
        #Make sure the input type is correct
        assert type(wire) == Wire.Wire
        #Set the values
        wire.setTo(self)
        self.fanIn.append(wire)
    
    def removeFanInWire(self,wire):
        '''
        removes a wire from fanIn. Not to be used by itself. Only to be called
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
        removes a wire from fanOut. Not to be used by itself. Only to be called
        by wires.
        '''
        assert type(wire) == Wire.Wire
        self.fanOut.remove(wire)
        wire.setFrom(None)
        
    def setFanOut(self,fanOut):
        '''
        Allows you to reset the fanOut to either an empty list or a list of wires
        '''
        #Inputs cannot have a fanOut
        if self.gateType == "Output":
            print "An Output gate cannot have a fanOut"
            return
        #Make sure the input type is correct
        assert type(fanOut) == list
        for wire in fanOut:
            assert type(wire) == Wire.Wire
        #Reset the To for the wires that were already in the fanOut
        for wire in self.fanOut:
            wire.setFrom(None)
        #set the values
        self.fanOut = []
        while len(fanOut)!=0:
            wire.setFrom(self)
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
        
    def setmB(self,mB):
        '''
        To Change the mB value
        '''
        self.mB = float(mB)
        
    def setpB(self,pB):
        '''
        To Change the mB value
        '''
        self.pB = float(pB)
        
    def setmRNAHalfLife(self,m_halfLife):
        '''
        set the mRNA half life and update the mRNA degradation rate
        '''
        if m_halfLife == None or (type(m_halfLife)!=int and type(m_halfLife)!=float):
            print "This is not a valid half life"
            return
        self.m_halfLife = float(m_halfLife)
        self.am = np.log(2)/self.m_halfLife
        
    def setProteinHalfLife(self,p_halfLife):
        '''
        set the protein half life and update the protein degradation rate
        '''
        if p_halfLife == None or (type(p_halfLife)!=int and type(p_halfLife)!=float):
            print "This is not a valid half life"
            return
        self.p_halfLife = float(p_halfLife)
        self.ap = np.log(2)/self.p_halfLife
        
    def setmRNADegRate(self,am):
        '''
        set the mRNA degradation rate and update the mRNA half life
        '''
        if am == None or (type(am)!=int and type(am)!=float):
            print "This is not a valid protein degradation rate"
            return
        self.am = float(am)
        self.m_halfLife = np.log(2)/self.am
        
    def setProteinDegRate(self,ap):
        '''
        set the protein degradation rate and update the protein half life
        '''
        if ap == None or (type(ap)!=int and type(ap)!=float):
            print "This is not a valid protein degradation rate"
            return
        self.ap = float(ap)
        self.p_halfLife = np.log(2)/self.ap
    
    def setInputType(self, inputType):
        '''
        Set the input to an equation from inputs
        '''
        assert self.gateType=='Input' and type(inputType)==list
        self.availableInputTypes = ['inputs.sinInput', 'inputs.linInput','inputs.sawInput','inputs.squInput', 'inputs.stepFunction', 'inputs.stepInput']
        self.numInputs = [4,2,3,4,3,4]
        assert inputType[0] in self.availableInputTypes
        #Make sure we have the correct number of inputs
        i = self.availableInputTypes.index(inputType[0])
        assert len(inputType) == self.numInputs[i]+1
        
        for i in range(1,len(inputType)):
            inputType[i] = float(inputType[i])
        self.inputType = inputType
        
    def setInitialmRNA(self, Mi):
        assert self.gateType!='Input'
        self.Mi = Mi
        
    def setInitialProtein(self, Pi):
        assert self.gateType!='Input'
        self.Pi = Pi
        
    def setDist(self, d):
        self.dist = int(d)
        
    def setVisited(self, visited):
        assert type(visited)==bool
        self.visited = visited

    def setExpectedProtein(self,tv):
        assert type(tv) == list or tv == None
        self.expectedProtein = tv
        if tv!=None:
            if self.gateType == "Input":
                self.setExpectedPromoter(tv)
            elif self.gateType == "Repressor":
                self.setExpectedPromoter(self.invertTruth(tv))
            elif self.gateType == "Output":
                self.setExpectedPromoter(None)
        elif tv==None:
            self.setExpectedPromoter(tv)
            
        
    def setExpectedPromoter(self,tv):
        assert type(tv) == list or tv == None
        self.expectedPromoter = tv
        
    def invertTruth(self, tv):
        newTruthValue = []
        for i in tv:
            if i == "0":
                newTruthValue.append("1")
            elif i == "1":
                newTruthValue.append("0")
        return newTruthValue
        
    def getAllProperties(self):
        '''
        Returns a dictionary of all the properties of the gate
        '''
        allProperties = {}
        allProperties['name'] = self.name
        allProperties['fanIn'] = self.fanIn
        allProperties['fanOut'] = self.fanOut
        allProperties['Km'] = self.Km
        allProperties['n'] = self.n
        allProperties['gateType'] = self.gateType
        allProperties['mB'] = self.mB
        allProperties['pB'] = self.pB
        allProperties['m_halfLife'] = self.m_halfLife
        allProperties['p_halfLife'] = self.p_halfLife
        allProperties['am'] = self.am
        allProperties['ap'] = self.ap
        allProperties['DIST'] = self.dist
        allProperties['EXPECTED_PROTEIN'] = self.expectedProtein
        allProperties['EXPECTED_PROMOTER'] = self.expectedPromoter
        if self.gateType=='Input':
            allProperties['InputType'] = self.inputType
        else:
            allProperties['Mi'] = self.Mi
            allProperties['Pi'] = self.Pi
        return allProperties
    
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
    def getmB(self):
        return self.mB
    def getpB(self):
        return self.pb
    def getm_halfLife(self):
        return self.m_halfLife
    def getp_halfLife(self):
        return self.p_halfLife
    def getam(self):
        return self.am
    def getap(self):
        return self.ap
    def getInputType(self):
        return self.inputType
    def getMi(self):
        return self.Mi
    def getPi(self):
        return self.Pi
    def getDist(self):
        return self.dist
    def wasVisited(self):
        return self.visited
    def getExpectedProtein(self):
        return self.expectedProtein
    def getExpectedPromoter(self):
        return self.expectedPromoter

    def formatJson(self):
        gateDescription = {}
        gateDescription['NAME'] = self.name
        gateDescription['DIST'] = self.dist
        gateDescription['EXPECTED_PROTEIN'] = self.expectedProtein
        gateDescription['EXPECTED_PROMOTER'] = self.expectedPromoter
        if self.gateType=='Input':
            gateDescription['TYPE'] = 'INPUT'
            gateDescription['INPUT'] = self.inputType
            
        elif self.gateType=='Repressor':
            if len(self.fanIn)==1:
                gateDescription['TYPE'] = 'NOT'
                gateDescription['INPUT'] = self.fanIn[0].getFrom().getName()
                gateDescription['Mi'] = self.Mi
                gateDescription['Pi'] = self.Pi
                gateDescription['m_HALFLIFE'] = self.m_halfLife
                gateDescription['p_HALFLIFE'] = self.p_halfLife
                gateDescription['am'] = 'NONE'
                gateDescription['mB'] = self.fanIn[0].getFrom().getmB()
                gateDescription['n'] = self.fanIn[0].getFrom().getn()
                gateDescription['Km'] = self.fanIn[0].getFrom().getKm()
                gateDescription['ap'] = 'NONE'
                gateDescription['pB'] = self.pB
                if self.fanIn[0].getFrom().getGateType()=="Input":
                    gateDescription['INPUT_EFFECT'] = "ACTIVATE"
                elif self.fanIn[0].getFrom().getGateType()=="Repressor":
                    gateDescription['INPUT_EFFECT'] = "REPRESS"
                
            elif len(self.fanIn)==2:
                gateDescription['TYPE'] = 'NOR'
                gateDescription['INPUT1'] = self.fanIn[0].getFrom().getName()
                gateDescription['INPUT2'] = self.fanIn[1].getFrom().getName()
                gateDescription['Mi'] = self.Mi
                gateDescription['Pi'] = self.Pi
                gateDescription['m_HALFLIFE'] = self.m_halfLife
                gateDescription['p_HALFLIFE'] = self.p_halfLife
                gateDescription['am'] = 'NONE'
                gateDescription['mB1'] = self.fanIn[0].getFrom().getmB()
                gateDescription['n1'] = self.fanIn[0].getFrom().getn()
                gateDescription['Km1'] = self.fanIn[0].getFrom().getKm()
                gateDescription['mB2'] = self.fanIn[1].getFrom().getmB()
                gateDescription['n2'] = self.fanIn[1].getFrom().getn()
                gateDescription['Km2'] = self.fanIn[1].getFrom().getKm()
                gateDescription['ap'] = 'NONE'
                gateDescription['pB'] = self.pB   
                if self.fanIn[0].getFrom().getGateType()=="Input":
                    gateDescription['INPUT1_EFFECT'] = "ACTIVATE"
                elif self.fanIn[0].getFrom().getGateType()=="Repressor":
                    gateDescription['INPUT1_EFFECT'] = "REPRESS"
                if self.fanIn[1].getFrom().getGateType()=="Input":
                    gateDescription['INPUT2_EFFECT'] = "ACTIVATE"
                elif self.fanIn[1].getFrom().getGateType()=="Repressor":
                    gateDescription['INPUT2_EFFECT'] = "REPRESS"                
                
        elif self.gateType=='Output':
            if len(self.fanIn)==1:
                gateDescription['TYPE'] = 'BUFFER'
                gateDescription['INPUT'] = self.fanIn[0].getFrom().getName()
                gateDescription['Mi'] = self.Mi
                gateDescription['Pi'] = self.Pi
                gateDescription['m_HALFLIFE'] = self.m_halfLife
                gateDescription['p_HALFLIFE'] = self.p_halfLife
                gateDescription['am'] = 'NONE'
                gateDescription['mB'] = self.fanIn[0].getFrom().getmB()
                gateDescription['n'] = self.fanIn[0].getFrom().getn()
                gateDescription['Km'] = self.fanIn[0].getFrom().getKm()
                gateDescription['ap'] = 'NONE'
                gateDescription['pB'] = self.pB
                if self.fanIn[0].getFrom().getGateType()=="Input":
                    gateDescription['INPUT_EFFECT'] = "ACTIVATE"
                elif self.fanIn[0].getFrom().getGateType()=="Repressor":
                    gateDescription['INPUT_EFFECT'] = "REPRESS"
                
                
            elif len(self.fanIn)==2:
                gateDescription['TYPE'] = 'OR'
                gateDescription['INPUT1'] = self.fanIn[0].getFrom().getName()
                gateDescription['INPUT2'] = self.fanIn[1].getFrom().getName()
                gateDescription['Mi'] = self.Mi
                gateDescription['Pi'] = self.Pi
                gateDescription['m_HALFLIFE'] = self.m_halfLife
                gateDescription['p_HALFLIFE'] = self.p_halfLife
                gateDescription['am'] = 'NONE'
                gateDescription['mB1'] = self.fanIn[0].getFrom().getmB()
                gateDescription['n1'] = self.fanIn[0].getFrom().getn()
                gateDescription['Km1'] = self.fanIn[0].getFrom().getKm()
                gateDescription['mB2'] = self.fanIn[1].getFrom().getmB()
                gateDescription['n2'] = self.fanIn[1].getFrom().getn()
                gateDescription['Km2'] = self.fanIn[1].getFrom().getKm()
                gateDescription['ap'] = 'NONE'
                gateDescription['pB'] = self.pB
                if self.fanIn[0].getFrom().getGateType()=="Input":
                    gateDescription['INPUT1_EFFECT'] = "ACTIVATE"
                elif self.fanIn[0].getFrom().getGateType()=="Repressor":
                    gateDescription['INPUT1_EFFECT'] = "REPRESS"
                if self.fanIn[1].getFrom().getGateType()=="Input":
                    gateDescription['INPUT2_EFFECT'] = "ACTIVATE"
                elif self.fanIn[1].getFrom().getGateType()=="Repressor":
                    gateDescription['INPUT2_EFFECT'] = "REPRESS"
        return gateDescription
        
    def __str__(self):
        return str(self.getAllProperties())
  