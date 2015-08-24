# -*- coding: utf-8 -*-
import LogicOperations

class Gate(object):
    def __init__(self,name,gateType,state="0"):
        assert type(name) == str
        assert gateType in ['INPUT','NOT','NOR','OR','AND','NAND','XNOR','XOR','BUF']
        self.name = name
        self.gateType = gateType
        self.state=str(state)
        self.fanOut = []
        self.fanIn = []
        self.dist = 0
        self.visited = False
        self.seenInput = False
        if gateType=="INPUT":
            self.seenInput= True

        
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
    def setCurrentState(self,state):
        state = str(state)
        assert state=="1" or state=="0"
        self.state = state
        
    def addFanOut(self,gate):
        """
        This is only used internally. Not intended to be called 
        by itself.
        """
        assert type(gate)==Gate
        self.fanOut.append(gate)
    def addFanIn(self,gate):
        """
        Adds a gate to the fanIn. Also adds this gate to the fanOut
        of the other gate.
        """
        assert type(gate)==Gate
        if self.gateType == "INPUT":
            print "Inputs cannot have a fanIn."
            return
        self.fanIn.append(gate)
        gate.addFanOut(self)
    def removeFanOut(self,gate):
        """
        This is only used internally. Not intended to be called 
        by itself.
        """
        assert type(gate)==Gate
        self.fanOut.remove(gate)
    def removeFanIn(self,gate):
        """
        Removes a gate from the fanIn of this one and removes this
        from the fanOut of the other gate.
        """
        assert type(gate)==Gate
        self.fanIn.remove(gate)
        gate.removeFanOut(self)
    def setVisited(self,status):
        self.visited = bool(status)
    def setSeenInput(self,status):
        self.seenInput = bool(status)
    def setDist(self,dist):
        self.dist = int(dist)
        
    #Getter functions
    def getName(self):
        return self.name
    def getGateType(self):
        return self.gateType
    def getFanIn(self):
        return self.fanIn
    def getFanOut(self):
        return self.fanOut
    def getCurrentState(self):
        return self.state
    def getVisited(self):
        return self.visited
    def getSeenInput(self):
        return self.seenInput
    def getDist(self):
        return self.dist
    def getNextState(self):
        operations = ['INPUT','NOT','NOR','OR','AND','NAND','XNOR','XOR','BUF']
        functions = [None,LogicOperations.NOT,LogicOperations.NOR,LogicOperations.OR,\
        LogicOperations.AND,LogicOperations.NAND,LogicOperations.XNOR,LogicOperations.XOR,LogicOperations.BUF]
        func = functions[operations.index(self.gateType)]
        tvList = []
        for gate in self.fanIn:
            tvList.append(gate.getCurrentState())
        if len(tvList)==0:
            return self.state
        return func(tvList)
    
    def __str__(self):
        """
        prints the gate in the way it would appear in a netlist
        """
        gateStr = self.gateType + "(" + self.name + ","
        for gate in self.fanIn:
            gateStr = gateStr + gate.getName() + ","
        #Remove the ending comma"
        gateStr = gateStr[:-1] + ")"
        return gateStr

class Graph(object):
    def __init__(self):
        self.isSequentialprop = False
        self.allGates = []
        self.inputGates = []
        self.otherGates = []
        self.history = []
        self.recentHistory = []
        
    def update(self):
        newState = []
        for gate in self.allGates:
            newState.append(gate.getNextState())
        for i in range(len(newState)):
            self.allGates[i].setCurrentState(newState[i])
        
        self.addToHistory()
        self.addToRecentHistory()
            
    def update2(self):
        for gate in self.allGates:
            gate.setCurrentState(gate.getNextState())        
        self.addToHistory()
        self.addToRecentHistory()

    def addGate(self,gate):
        assert type(gate)==Gate
        self.allGates.append(gate)
        if gate.getGateType()=="INPUT":
            self.inputGates.append(gate)
        else:
            self.otherGates.append(gate)
        #Reset history everytime you add or remove gate
        self.history = [self.getState()]
        self.recentHistory = [self.getState()]
        self.checkSequential()
        self.otherGates.sort(self.gateCompare)
        
    def removeGate(self,gate):
        assert type(gate)==Gate
        if gate in self.allGates:
            self.allGates.remove(gate)
            if gate.getGateType()=="INPUT":
                self.inputGates.remove(gate)
            else:
                self.otherGates.remove(gate)
            #Reset history everytime you add or remove gate
            self.history = [self.getState()]
            self.recentHistory = [self.getState()]
        else:
            print "This gate is not in the graph"
    def getState(self):
        state = []
        for gate in self.allGates:
            state.append(gate.getCurrentState())
        return state
        
    def addToHistory(self):
        self.history.append(self.getState())

    def addToRecentHistory(self):
        self.recentHistory.append(self.getState())
        
    def alreadySeenState(self):
        """
        Returns a tuple where the first element tells whether or not
        we have seen the current state and the second element tells
        whether it is a repeat of the last state meaning it is stable
        """
        if self.recentHistory.count(self.getState())>1:
            if self.getState() == self.recentHistory[-2]:
                return (True,True)
            return (True,False)
        return (False,False)
        
    def printNames(self):
        names = ""
        for gate in self.allGates:
            names = names + gate.getName() + "\t"
        print names
        
    def printStates(self):
        states = ""
        for gate in self.allGates:
            states = states + gate.getCurrentState() + "\t"
            
    def runUntilStableOrRepeat(self,cond):
        """
        if cond==1, it will print full history
        if cond==0, it will print recent history
        otherwise it will not print any history
        """

        while not(self.alreadySeenState()[0]):
            #self.update()
            self.update2()
        if cond==1:
            self.printHistory()
        elif cond==0:
            self.printRecentHistory()
        if self.alreadySeenState()[1]:
            print "stable"
        else:
            print "unstable"
            
    def printHistory(self):
        self.printNames()
        for state in self.history:
            print '\t'.join(state)
            
    def printRecentHistory(self):
        self.printNames()
        for state in self.recentHistory:
            print '\t'.join(state)
            
    def getGate(self,name):
        #Given a name return the gate
        for gate in self.allGates:
            if gate.getName()==name:
                return gate
        print "No gate found by this name."
        return
        
    def setGateStates(self,states):
        assert type(states)==list and len(states)==len(self.allGates)
        for i in range(len(self.allGates)):
            self.allGates[i].setCurrentState(states[i])
        #self.addToHistory()
        self.addToHistory()
        self.recentHistory = [self.getState()]
    
    def checkSequential(self):
        self.isSequentialprop = False
        for gate in self.allGates:
            gate.setVisited(False)
            gate.setDist(0)
            if gate.getGateType()!="INPUT":
                gate.setSeenInput(False)
        for gate in self.allGates:
            self.recursivelyFindInputs(gate)
    
    def recursivelyFindInputs(self,gate):
        assert type(gate)==Gate
        if gate.getSeenInput():
            return True
        if gate.getVisited():
            self.isSequentialprop = True
            return False
        gate.setVisited(True)
        components = []
        toReturn = True
        seenInput = True
        for fanInGate in gate.getFanIn():
            
            if not(self.recursivelyFindInputs(fanInGate)):
                toReturn = False
                seenInput = False
            components.append(fanInGate.getDist())
        gate.setDist(max(components)+1)
        gate.setSeenInput(seenInput)
        return toReturn
    
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
    
    
    def isSequential(self):
        self.checkSequential()
        return self.isSequentialprop
            
    def __str__(self):
        fullString = ""
        for gate in self.allGates:
            fullString = fullString + str(gate) + " state= " + gate.getCurrentState() + "\n"
        return fullString
        
def test():
    gate1 = Gate("R","INPUT",0)
    gate2 = Gate("S","INPUT",0)
    gate3 = Gate("g3","NOR",0)
    gate4 = Gate("g4","NOR",0)
    
    gate3.addFanIn(gate1)
    gate3.addFanIn(gate4) 
    gate4.addFanIn(gate2)
    gate4.addFanIn(gate3)
    
    graph = Graph()
    graph.addGate(gate1)
    graph.addGate(gate2)
    graph.addGate(gate3)
    graph.addGate(gate4)
    
    print graph
    print graph.isSequential()
#    print graph
#    print "break"
#    graph.printHistory()
#    print "break"
#    graph.setGateStates([1,0,0,0])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([0,1,0,0])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([0,1,1,0])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([0,1,0,1])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
    
def test2():
    gate1 = Gate("g1","NOT",0)
    gate2 = Gate("g2","NOT",0)
    gate3 = Gate("g3","NOT",0)
    
    gate1.addFanIn(gate2)
    gate2.addFanIn(gate3) 
    gate3.addFanIn(gate1)
    
    graph = Graph()
    graph.addGate(gate1)
    graph.addGate(gate2)
    graph.addGate(gate3)
    
    print graph.isSequential()
    print graph
    print "break"
#    graph.printHistory()
#    print "break"
#    graph.setGateStates([0,0,0])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([1,0,0])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([0,1,0])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([0,0,1])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([0,1,1])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
#    graph.setGateStates([1,0,1])
#    graph.runUntilStableOrRepeat(0)
#    print "break"
    graph.setGateStates([1,1,0])
    graph.runUntilStableOrRepeat(0)
    print "break"
#    graph.setGateStates([1,1,1])
#    graph.runUntilStableOrRepeat(0)
#    print "break"