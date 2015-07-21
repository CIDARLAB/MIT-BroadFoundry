import random
import time

#Helper functions for counting gates and checking truth value
def getTruthValue(circuit):
    circuit = circuit.replace("0","00000000")
    circuit = circuit.replace("1","11111111")
    circuit = circuit.replace("a'","11110000")
    circuit = circuit.replace("b'","11001100")
    circuit = circuit.replace("c'","10101010")
    circuit = circuit.replace("a","00001111")
    circuit = circuit.replace("b","00110011")
    circuit = circuit.replace("c","01010101")
    while len(circuit) != 8:
        indeces = findSmallestParentheses(circuit)
        group = circuit[indeces[0]:indeces[1]+1]
        shouldInvert = False
        try:
            char = circuit[indeces[1]+1]
            if char=="'":
                shouldInvert = True
        except:
            pass
        if group[9] == "&":
            answer = andIt(group)
        elif group[9] == ".":
            answer = nor(group)
        if shouldInvert:
            group = circuit[indeces[0]:indeces[1]+2]
            answer = invert(answer)
            circuit = circuit[0:indeces[0]] + answer + circuit[indeces[1]+2:]
        elif not shouldInvert:
            circuit = circuit[0:indeces[0]] + answer + circuit[indeces[1]+1:]
##        circuit = circuit.replace(group,answer)
    return circuit

def gateCounter(circuit):
    indeces = findLeftParentheses(circuit)
    tokens = []
    while indeces != None:
        if tokens.count(circuit[indeces[0]:indeces[1]+1]) == 0:
            tokens.append(circuit[indeces[0]:indeces[1]+1])
        circuit = circuit[0:indeces[0]] + circuit[indeces[0]+1:indeces[1]]+circuit[indeces[1]+1:len(circuit)]
        
        indeces = findLeftParentheses(circuit)
    return len(tokens)

#Other helper functions
def findSmallestParentheses(circuit):
    #Assuming valid pairing of parentheses
    L = len(circuit)
    for i in xrange(L):
        if circuit[i]=="(":
            startIndex = i
        if circuit[i]==")":
            endIndex = i
            return (startIndex,endIndex)
    return None
def nor(s):
    a = s[1:9]
    b = s[10:18]
    answer = ""
    for i in xrange(8):
        if (a[i]=="1" or b[i]=="1"):
            answer += "0"
        elif (a[i]=="0"  and b[i]=="0"):
            answer += "1"
    return answer
def andIt(s):
    a = s[1:9]
    b = s[10:18]
    answer = ""
    for i in xrange(8):
        if (a[i]=="1" and b[i]=="1"):
            answer += "1"
        else:
            answer += "0"
    return answer
def invert(s):
    answer = ""
    for i in xrange(8):
        if (s[i]=="1"):
            answer += "0"
        else:
            answer += "1"
    return answer
def findLeftParentheses(circuit):
    L = len(circuit)
    counter = 0
    indeces = []
    for i in xrange(L):
        if circuit[i] == "(":
            counter += 1
            if len(indeces)==0:
                indeces.append(i)
        elif circuit[i] == ")":
            counter -= 1
            if counter == 0:
                indeces.append(i)
                return indeces
    return None

def findLeftSingleQuotes(s):
    L = len(s)
    counter = 0
    indeces = []
    for i in xrange(L):
        if s[i] == "'":
            if len(indeces)==0:
                indeces.append(i)
            elif len(indeces)==1:
                indeces.append(i)
                return indeces
    return None
def findLeftDoubleQuotes(s):
    L = len(s)
    counter = 0
    indeces = []
    for i in xrange(L):
        if s[i] == '"':
            if len(indeces)==0:
                indeces.append(i)
            elif len(indeces)==1:
                indeces.append(i)
                return indeces
    return None

def sortAndPrint(foundCircuits):
    truths = foundCircuits.keys()
    truths.sort()
    for truth in truths:
        print truth," : ",foundCircuits[truth]
        
def getFromFile(fileName):
    foundCircuits = {}
    myFile = open(fileName,'r')
    for line in myFile:
        #get the truth value
        key = line[0:8]
        #prepare the list of circuits
        value = []
        #get the rest of the line removing the \n
        unformatedValue = line[9:-1]
        #remove the surrounding brackets
        unformatedValue = unformatedValue[1:-1]
        
        #For the circuits with ' in them
        doubleQuoteIndeces = findLeftDoubleQuotes(unformatedValue)
        while doubleQuoteIndeces != None:
            #we dont want to include the quotation marks themselves
            value.append(unformatedValue[doubleQuoteIndeces[0]+1:doubleQuoteIndeces[1]])
            unformatedValue = unformatedValue[0:doubleQuoteIndeces[0]]+unformatedValue[doubleQuoteIndeces[1]+1:len(unformatedValue)]
            doubleQuoteIndeces = findLeftDoubleQuotes(unformatedValue)
            
        #For the circuits without ' in them
        singleQuoteIndeces = findLeftSingleQuotes(unformatedValue)
        while singleQuoteIndeces != None:
            #we dont want to include the quotation marks themselves
            value.append(unformatedValue[singleQuoteIndeces[0]+1:singleQuoteIndeces[1]])
            unformatedValue = unformatedValue[singleQuoteIndeces[1]+1:len(unformatedValue)]
            singleQuoteIndeces = findLeftSingleQuotes(unformatedValue)
        #add it to the dictionary
        foundCircuits[line[0:8]] = value
    myFile.close()
    print len(foundCircuits)
    return foundCircuits

def writeToFile(foundCircuits, fileName):
    myFile = open(fileName,'w')
    keys = foundCircuits.keys()
    keys.sort()
    for key in keys:
        s = key + " " + str(foundCircuits[key]) + "\n"
        myFile.write(s)
    myFile.close()
        
def checkForUse(gateCount,notUsed,fileName):
    allFoundCircuits = minCircuitFinder2(gateCount,notUsed)[1]
    truthCircuits = getFromFile(fileName)
    unused = []
    used = []
    printed = []
    print "start"
    for circuit in allFoundCircuits[gateCount]:
        for listOfCircuits in truthCircuits.values():
            for usedCircuit in listOfCircuits:
                if usedCircuit.count(circuit)>0:
                    if used.count(circuit)==0:
                        used.append(circuit)
                        #Don't know how to break just the inner two for loops so I only break one
                        break
    for circuit in allFoundCircuits[gateCount]:
        if used.count(circuit)==0:
            unused.append(circuit)
    print "Of the "+str(len(allFoundCircuits[gateCount]))+" total cicuits with "+str(gateCount)+" gates,"
    print str(len(unused))+" unused circuits with "+str(gateCount)+" gates"
    return unused#,used   

def minNumberOfCircuitsFor256():
    abc = [1,2,3,4,5,6,7,8]
    acb = [1,3,2,4,5,7,6,8]
    bac = [1,2,5,6,3,4,7,8]
    bca = [1,5,2,6,3,7,4,8]
    cab = [1,3,5,7,2,4,6,8]
    cba = [1,5,3,7,2,6,4,8]
    needed = []
    
    for i in range(256):
        x = "{0:b}".format(i)
        while len(x)<8:
            x = "0"+x
        needed.append(x)
    madeFrom = {}
    found = []
    for truth in needed:
        if truth not in found:
            madeFrom[truth]=[truth]
            acbtruth = ""
            bactruth = ""
            bcatruth = ""
            cabtruth = ""
            cbatruth = ""
            for i in range(8):
                acbtruth += truth[abc[acb[i]-1]-1]
                bactruth += truth[abc[bac[i]-1]-1]
                bcatruth += truth[abc[bca[i]-1]-1]
                cabtruth += truth[abc[cab[i]-1]-1]
                cbatruth += truth[abc[cba[i]-1]-1]
            
            madeFrom[truth].append(acbtruth)
            madeFrom[truth].append(bactruth)
            madeFrom[truth].append(bcatruth)
            madeFrom[truth].append(cabtruth)
            madeFrom[truth].append(cbatruth)

            found.append(truth)
            if acbtruth not in found:
                found.append(acbtruth)
            if bactruth not in found:
                found.append(bactruth)
            if bcatruth not in found:
                found.append(bcatruth)
            if cabtruth not in found:
                found.append(cabtruth)
            if cbatruth not in found:
                found.append(cbatruth)
            
            
    print str(len(madeFrom)) + " circuits needed."
    return madeFrom

def findIsomorphs(circuit):
    type1 = circuit
    type2 = circuit
    type3 = circuit
    type4 = circuit
    type5 = circuit
    
    type1.replace("b","x")
    type1.replace("c","b")
    type1.replace("x","c")

    type2.replace("b","x")
    type2.replace("a","b")
    type2.replace("x","a")

    type3.replace("a","x")
    type3.replace("c","a")
    type3.replace("b","c")
    type3.replace("x","b")

    type4.replace("a","x")
    type4.replace("b","a")
    type4.replace("c","b")
    type4.replace("x","c")

    type5.replace("a","x")
    type5.replace("c","a")
    type5.replace("x","c")
    isomorphs = [circuit, type1,type2,type3,type4,type5]
    for isomorph in isomorphs:
        isomorph.replace("(a.b)","(b.a)")
        isomorph.replace("(a.c)","(c.a)")
        isomorph.replace("(a.0)","(0.a)")
        isomorph.replace("(b.c)","(c.b)")
        isomorph.replace("(b.0)","(0.b)")
        isomorph.replace("(c.0)","(0.c)")
    return isomorphs

def toBinary(num):
    #converts a decimal number to an 8 bit binary string
    s = "{0:b}".format(num)
    while len(s)<8:
        s = "0"+s
    return s

def containsEqualSubCircuit(circuit):
    truthTable = getTruthValue(circuit)
    circuit = circuit.replace("0","00000000")
    circuit = circuit.replace("1","11111111")
    circuit = circuit.replace("a'","11110000")
    circuit = circuit.replace("b'","11001100")
    circuit = circuit.replace("c'","10101010")
    circuit = circuit.replace("a","00001111")
    circuit = circuit.replace("b","00110011")
    circuit = circuit.replace("c","01010101")
    while len(circuit) != 8:
        if circuit.find(truthTable)>0:
            return True
        indeces = findSmallestParentheses(circuit)
        group = circuit[indeces[0]:indeces[1]+1]
        shouldInvert = False
        try:
            char = circuit[indeces[1]+1]
            if char=="'":
                shouldInvert = True
        except:
            pass
        if group[9] == "&":
            answer = andIt(group)
        elif group[9] == ".":
            answer = nor(group)
        if shouldInvert:
            group = circuit[indeces[0]:indeces[1]+2]
            answer = invert(answer)
        circuit = circuit.replace(group,answer)
    return False


#CircuitFinders
def minCircuitFinder1(num):
    allLevels = {}
    currentLevel = ['a','b','c',"a'","b'","c'"]
    allLevels[0] = currentLevel
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = ['a']
    foundCircuits[getTruthValue('b')] = ['b']
    foundCircuits[getTruthValue('c')] = ['c']
    foundCircuits[getTruthValue("a'")] = ["a'"]
    foundCircuits[getTruthValue("b'")] = ["b'"]
    foundCircuits[getTruthValue("c'")] = ["c'"]
    foundCircuits[getTruthValue('0')] = ['0']
    foundCircuits[getTruthValue('1')] = ['1']
    reallyNotAllowed = [getTruthValue('0'),getTruthValue('1')]
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue("a'"),getTruthValue("b'"),getTruthValue("c'"),
                  getTruthValue("(a&b)"),getTruthValue("(a&c)"),getTruthValue("(b&c)"),
                  getTruthValue("(a&b)'"),getTruthValue("(a&c)'"),getTruthValue("(b&c)'"),
                  getTruthValue("(a&b')"),getTruthValue("(a&c')"),getTruthValue("(b&c')"),
                  getTruthValue("(a&b')'"),getTruthValue("(a&c')'"),getTruthValue("(b&c')'"),
                  getTruthValue("(a'&b)"),getTruthValue("(a'&c)"),getTruthValue("(b'&c)"),
                  getTruthValue("(a'&b)'"),getTruthValue("(a'&c)'"),getTruthValue("(b'&c)'"),
                  getTruthValue("(a'&b')"),getTruthValue("(a'&c')"),getTruthValue("(b'&c')"),
                  getTruthValue("(a'&b')'"),getTruthValue("(a'&c')'"),getTruthValue("(b'&c')'"),
                  ]
    startTime = time.time()
    levelGates=0
    while not(levelGates>=num+1 or len(foundCircuits)>=256):
        foundTruths = []
        currentLevel = allLevels[levelGates]
        for alphaCircuit in currentLevel:
            alphaTruthTable = getTruthValue(alphaCircuit)
            if not(alphaTruthTable in foundCircuits):
                foundCircuits[alphaTruthTable] = [alphaCircuit]
                foundTruths.append(getTruthValue(alphaCircuit))
            elif alphaTruthTable in foundCircuits:
                if alphaTruthTable in foundTruths:
                    foundCircuits[alphaTruthTable].append(alphaCircuit)
        if len(foundCircuits)>=256:
            return foundCircuits,allLevels
        if levelGates == num:
            endTime = time.time()
            totalTime = endTime-startTime
            print "It took "+str(totalTime)+" to finish."
            return foundCircuits,allLevels
        print "Found " + str(len(foundCircuits))+" truthValues found so far."
        for alphaCircuit in currentLevel:
            for comparingLevelGates in range(levelGates+1):
                comparingLevel = allLevels[comparingLevelGates]
                if comparingLevel==currentLevel:
                    for betaCircuitIndex in range(comparingLevel.index(alphaCircuit)+1,len(comparingLevel)):
                        betaCircuit = comparingLevel[betaCircuitIndex]
                        if betaCircuit == alphaCircuit:
                            break
                        newCircuit = "("+alphaCircuit+"&"+betaCircuit+")"
                        newCircuit2 = newCircuit + "'"
                        numGates = gateCounter(newCircuit)
                        if (numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit)
                            except KeyError:
                                allLevels[numGates] = [newCircuit]
                        if (numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit2)
                            except KeyError:
                                allLevels[numGates] = [newCircuit2]
                    break
                elif comparingLevel!=currentLevel:
                    for betaCircuit in comparingLevel:
                        newCircuit = "("+alphaCircuit+"&"+betaCircuit+")"
                        newCircuit2 = newCircuit + "'"
                        numGates = gateCounter(newCircuit)
                        if (numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit)
                            except KeyError:
                                allLevels[numGates] = [newCircuit]
                        if (numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit2)
                            except KeyError:
                                allLevels[numGates] = [newCircuit2]
        print "Done with "+str(levelGates)+" gate circuits"
        print str(levelGates+1)+" has "+str(len(allLevels[levelGates+1]))+" circuits."
        endTime = time.time()
        totalTime = endTime-startTime
        print "It took "+str(totalTime)+" to get this far."

        levelGates += 1
##        for key in allLevels:
##            print str(key)+" has "+str(len(allLevels[key]))+" circuits."
##        pause = raw_input("Done with "+str(levelGates))
    return foundCircuits, allLevels

def minCircuitFinder2(num,ignore):
    #ignore is a list of circuits to ignore
    allLevels = {}
    currentLevel = ['a','b','c',"a'","b'","c'"]
    allLevels[0] = currentLevel
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = ['a']
    foundCircuits[getTruthValue('b')] = ['b']
    foundCircuits[getTruthValue('c')] = ['c']
    foundCircuits[getTruthValue("a'")] = ["a'"]
    foundCircuits[getTruthValue("b'")] = ["b'"]
    foundCircuits[getTruthValue("c'")] = ["c'"]
    foundCircuits[getTruthValue('0')] = ['0']
    foundCircuits[getTruthValue('1')] = ['1']
    reallyNotAllowed = [getTruthValue('0'),getTruthValue('1')]
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue("a'"),getTruthValue("b'"),getTruthValue("c'"),
                  getTruthValue("(a&b)"),getTruthValue("(a&c)"),getTruthValue("(a&0)"),
                  getTruthValue("(b&c)"),getTruthValue("(b&0)"),getTruthValue("(c&0)"),
                  getTruthValue("(a&b)'"),getTruthValue("(a&c)'"),getTruthValue("(a&0)'"),
                  getTruthValue("(b&c)'"),getTruthValue("(b&0)'"),getTruthValue("(c&0)'")]
    startTime = time.time()
    levelGates=0
    while not(levelGates>=num+1 or len(foundCircuits)>=256):
        foundTruths = []
        currentLevel = allLevels[levelGates]
        for alphaCircuit in currentLevel:
            alphaTruthTable = getTruthValue(alphaCircuit)
            if not(alphaTruthTable in foundCircuits):
                foundCircuits[alphaTruthTable] = [alphaCircuit]
                foundTruths.append(getTruthValue(alphaCircuit))
            elif alphaTruthTable in foundCircuits:
                if getTruthValue(alphaCircuit) in foundTruths:
                    foundCircuits[alphaTruthTable].append(alphaCircuit)
        if len(foundCircuits)>=256:
            return foundCircuits,allLevels
        if levelGates == num:
            endTime = time.time()
            totalTime = endTime-startTime
            print "It took "+str(totalTime)+" to finish."
            return foundCircuits,allLevels
        print "Found " + str(len(foundCircuits))+" truthValues found so far."
        for alphaCircuit in currentLevel:
            for comparingLevelGates in range(levelGates+1):
                comparingLevel = allLevels[comparingLevelGates]
                if comparingLevel==currentLevel:
                    for betaCircuitIndex in range(comparingLevel.index(alphaCircuit)+1,len(comparingLevel)):
                        betaCircuit = comparingLevel[betaCircuitIndex]
                        if betaCircuit == alphaCircuit:
                            break
                        newCircuit = "("+alphaCircuit+"&"+betaCircuit+")"
                        newCircuit2 = newCircuit + "'"
                        numGates = gateCounter(newCircuit)
                        if ((newCircuit in ignore)
                            or numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit)
                            except KeyError:
                                allLevels[numGates] = [newCircuit]
                        if ((newCircuit in ignore)
                            or numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit2)
                            except KeyError:
                                allLevels[numGates] = [newCircuit2]
                    break
                elif comparingLevel!=currentLevel:
                    for betaCircuit in comparingLevel:
                        newCircuit = "("+alphaCircuit+"&"+betaCircuit+")"
                        newCircuit2 = newCircuit + "'"
                        numGates = gateCounter(newCircuit)
                        if ((newCircuit in ignore)
                            or numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit)
                            except KeyError:
                                allLevels[numGates] = [newCircuit]
                        if ((newCircuit in ignore)
                            or numGates>num
                            or (numGates>1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in notAllowed)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))
                            or (numGates==1
                                and(containsEqualSubCircuit(newCircuit)
                                    or (getTruthValue(newCircuit) in reallyNotAllowed)))):
                            pass
                        else:
                            try:
                                allLevels[numGates].append(newCircuit2)
                            except KeyError:
                                allLevels[numGates] = [newCircuit2]
        print "Done with "+str(levelGates)+" gate circuits"
        print str(levelGates+1)+" has "+str(len(allLevels[levelGates+1]))+" circuits."
        endTime = time.time()
        totalTime = endTime-startTime
        print "It took "+str(totalTime)+" to get this far."

        levelGates += 1
##        for key in allLevels:
##            print str(key)+" has "+str(len(allLevels[key]))+" circuits."
##        pause = raw_input("Done with "+str(levelGates))
    return foundCircuits, allLevels

def wrapWithUnused(startNumGates,endNumGates,directory):
    #make the first level
    unused = []
    gates = minCircuitFinder2(startNumGates,[])
    if directory[-1] != "/":
        directory += "/"
    startFileName = directory+"start"+str(startNumGates)+"gates.txt"
    writeToFile(gates[0],startFileName)
    for i in range(startNumGates+1):
        unused = unused + checkForUse(i,unused,startFileName)
    for currentGates in range(startNumGates+1,endNumGates+1):
        if currentGates < endNumGates:
            currentFileName = directory+str(currentGates)+"gates.txt"
            gates = minCircuitFinder2(currentGates,unused)
            writeToFile(gates[0],currentFileName)
            unused = unused + checkForUse(currentGates,unused,currentFileName)
        if currentGates == endNumGates:
            currentFileName = directory+str(currentGates)+"gates.txt"
            gates = minCircuitFinder2(currentGates,unused)
            writeToFile(gates[0],currentFileName)
            return gates
    return gates


def wrapWithUnusedFromStock(startNumGates,endNumGates,directory):
    #use the stock
    startTime = time.time()
    unused = []
    directory.replace("\\","/")
    if directory[-1] != "/":
        directory += "/"
    startFileName = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/GateMinimizationAND/Stock/"+str(startNumGates)+"gates.txt"
    for i in range(startNumGates+1):
        unused = unused + checkForUse(i,unused,startFileName)
    for currentGates in range(startNumGates+1,endNumGates+1):
        if currentGates < endNumGates:
            currentFileName = directory+str(currentGates)+"gates.txt"
            gates = minCircuitFinder2(currentGates,unused)
            writeToFile(gates[0],currentFileName)
            unused = unused + checkForUse(currentGates,unused,currentFileName)
        if currentGates == endNumGates:
            currentFileName = directory+str(currentGates)+"gates.txt"
            gates = minCircuitFinder2(currentGates,unused)
            writeToFile(gates[0],currentFileName)
            return gates
    endTime = time.time()
    overallTime = endTime-startTime
    print "Overall this took " + str(overallTime) + " seconds."
    return gates
