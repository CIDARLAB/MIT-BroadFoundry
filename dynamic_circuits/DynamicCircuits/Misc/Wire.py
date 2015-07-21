# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 22:31:03 2015

@author: Arinze
"""
import numpy as np
import Gate
import DAG

class Wire(object):
    
    def __init__(self,name):
        assert type(name) == str
        self.name = name
        self.fromGate = None
        self.toGate = None
        
    def getName(self):
        return self.name
        
    def getFrom(self):
        return self.fromGate
        
    def getTo(self):
        return self.toGate

    def setName(self,name):
        assert type(name)==name
        self.name = name

    def setFrom(self, fromGate):
        '''
        Sets the gate at the end of the wire. This should not be used by itself.
        This is only to be used by the Gate object when we are setting the fanOut.
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
            
    def __eq__(self, wire):
        if type(wire) != Wire:
            return False
        if (self.toGate == wire.getTo()) and (self.fromGate == wire.getFrom()) and (self.name == wire.getName()):
            return True
        return False
    def __ne(self,wire):
        return not(self.__eq__(wire))
        
    def __str__(self):
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