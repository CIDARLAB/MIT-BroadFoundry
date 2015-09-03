# -*- coding: utf-8 -*-
"Real"
"""
Created on Fri Jul 17 22:31:03 2015

@author: Arinze
Edited by Alex on 7/30/2015
"""

import Gate

class Wire(object):
    
    def __init__(self,name):
        assert type(name) == str
        self.name = name
        self.fromGate = None
        self.toGate = None
        
    #getters
    def getName(self):
        return self.name
        
    def getFrom(self):
        return self.fromGate
        
    def getTo(self):
        return self.toGate
        
    #setters
    def setName(self,name):
        assert type(name)== str
        self.name = name

    def setFrom(self, fromGate):
        '''
        Sets the gate at the end of the wire. This should not be used by itself.
        This is only to be used by the Gate object when setting the fanOut.
        '''
        assert type(fromGate)== Gate.Gate or fromGate == None
        if self.fromGate != None:
            try:
                self.fromGate.removeFanOutWire(self)
            except ValueError:
                pass
        #Since this will only be called by the gate, the gate itself will change
        #its own fanOut
        self.fromGate = fromGate
        
    def setTo(self, toGate):
        '''
        Sets the gate at the end of the wire. This should not be used by itself.
        This is only to be used by the Gate object when we are setting the fanIn.
        '''
        assert type(toGate)== Gate.Gate or toGate == None
        if self.toGate != None:
            try:
                self.toGate.removeFanInWire(self)
            except ValueError:
                pass
        #Since this will only be called by the gate, the gate itself will change
        #its own fanIn
        self.toGate = toGate
            
        
    def __repr__(self,wire):
        return self.name
    
    def __str__(self):
        '''
        String format for wire
        example: 'g1--w1->g2'
        '''
        result = ''
        try:
            result = result + self.fromGate.getName()
        except AttributeError:
            result = result + "None"
        result = result + '--' + self.name + '->'
        try:
            result = result + self.toGate.getName()
        except AttributeError:
            result = result + "None"
        return result