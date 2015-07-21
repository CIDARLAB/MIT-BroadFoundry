import random
import time

def levelMaker1(numLevels):
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    currentLevel.append("("+alphaCircuit+"."+betaCircuit+")")
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                currentLevel.append("("+alphaCircuit+"."+betaCircuit+")")
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return allLevels
            
    
    
def getTruthValue(circuit):
    circuit = circuit.replace("0","00000000")
    circuit = circuit.replace("1","11111111")
    circuit = circuit.replace("a","00001111")
    circuit = circuit.replace("b","00110011")
    circuit = circuit.replace("c","01010101")
    while len(circuit) != 8:
        indeces = findSmallestParentheses(circuit)
        toNor = circuit[indeces[0]:indeces[1]+1]
        answer = nor(toNor)
        circuit = circuit.replace(toNor,answer)
    return circuit
    

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
    
def gateCounter(circuit):
    indeces = findLeftParentheses(circuit)
    tokens = []
    while indeces != None:
        if tokens.count(circuit[indeces[0]:indeces[1]+1]) == 0:
            tokens.append(circuit[indeces[0]:indeces[1]+1])
        circuit = circuit[0:indeces[0]] + circuit[indeces[0]+1:indeces[1]]+circuit[indeces[1]+1:len(circuit)]
        
        indeces = findLeftParentheses(circuit)
    return len(tokens)

def levelMaker2(numLevels):
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    numGates = gateCounter(newCircuit)

                    if (i==1 or notAllowed.count(truthTable)==0):
                        if i==1:
                            notAllowed.append(truthTable)
                        if (foundCircuits.has_key(truthTable)):
                            if numGates< gateCounter(foundCircuits[truthTable]):
                                foundCircuits[truthTable] = newCircuit
                        if (not(foundCircuits.has_key(truthTable))):
                            foundCircuits[truthTable] = newCircuit
                        currentLevel.append(newCircuit)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                numGates = gateCounter(newCircuit)

                if (i==1 or notAllowed.count(truthTable)==0):
                    if i==1:
                        notAllowed.append(truthTable)
                    if (foundCircuits.has_key(truthTable)):
                        if numGates< gateCounter(foundCircuits[truthTable]):
                            foundCircuits[truthTable] = newCircuit
                    if (not(foundCircuits.has_key(truthTable))):
                        foundCircuits[truthTable] = newCircuit
                    currentLevel.append(newCircuit)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker3(numLevels):
    #same as levelMaker2 but stops as soon as dictionary reaches 256
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    numGates = gateCounter(newCircuit)

                    if (i==1 or notAllowed.count(truthTable)==0):
                        if i==1:
                            notAllowed.append(truthTable)
                        if (foundCircuits.has_key(truthTable)):
                            if numGates< gateCounter(foundCircuits[truthTable]):
                                foundCircuits[truthTable] = newCircuit
                                if (len(foundCircuits)==256):
                                    return foundCircuits
                        if (not(foundCircuits.has_key(truthTable))):
                            foundCircuits[truthTable] = newCircuit
                            if (len(foundCircuits)==256):
                                return foundCircuits
                        currentLevel.append(newCircuit)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                numGates = gateCounter(newCircuit)

                if (i==1 or notAllowed.count(truthTable)==0):
                    if i==1:
                        notAllowed.append(truthTable)
                    if (foundCircuits.has_key(truthTable)):
                        if numGates< gateCounter(foundCircuits[truthTable]):
                            foundCircuits[truthTable] = newCircuit
                            if (len(foundCircuits)==256):
                                return foundCircuits
                    if (not(foundCircuits.has_key(truthTable))):
                        foundCircuits[truthTable] = newCircuit
                        if (len(foundCircuits)==256):
                            return foundCircuits
                    currentLevel.append(newCircuit)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker4(numLevels):
    #the same as levelMaker2 but does not create a dictionary as it goes
    #prints the size of the current level being built with a 1/10000 chance
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    chance = 0.0001
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    if (i==1 or notAllowed.count(truthTable)==0):
                        if i==1:
                            notAllowed.append(truthTable)
                        
                        currentLevel.append(newCircuit)
                        if random.random()<chance:
                            print len(currentLevel)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                if (i==1 or notAllowed.count(truthTable)==0):
                    if i==1:
                        notAllowed.append(truthTable)
                    currentLevel.append(newCircuit)
                    if random.random()<chance:
                        print len(currentLevel)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker5(numLevels):
    #same as levelMaker2 but stops as soon as dictionary reaches 256
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    chance = 0.0001
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    numGates = gateCounter(newCircuit)

                    if (i==1 or notAllowed.count(truthTable)==0):
                        if i==1:
                            notAllowed.append(truthTable)
                        if (foundCircuits.has_key(truthTable)):
                            if numGates< gateCounter(foundCircuits[truthTable]):
                                foundCircuits[truthTable] = newCircuit
                                if (len(foundCircuits)==256):
                                    return foundCircuits
                        if (not(foundCircuits.has_key(truthTable))):
                            foundCircuits[truthTable] = newCircuit
                            if (len(foundCircuits)==256):
                                return foundCircuits
                        currentLevel.append(newCircuit)
                        if random.random()<chance:
                            print "number of truthTables found: " + str(len(foundCircuits))
                            print len(currentLevel)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                numGates = gateCounter(newCircuit)

                if (i==1 or notAllowed.count(truthTable)==0):
                    if i==1:
                        notAllowed.append(truthTable)
                    if (foundCircuits.has_key(truthTable)):
                        if numGates< gateCounter(foundCircuits[truthTable]):
                            foundCircuits[truthTable] = newCircuit
                            if (len(foundCircuits)==256):
                                return foundCircuits
                    if (not(foundCircuits.has_key(truthTable))):
                        foundCircuits[truthTable] = newCircuit
                        if (len(foundCircuits)==256):
                            return foundCircuits
                    currentLevel.append(newCircuit)
                    if random.random()<chance:
                        print "number of truthTables found: " + str(len(foundCircuits))
                        print len(currentLevel)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker6(numLevels):
    #same as two but does not save the values of level 6
    #took >14hrs to get to 254/256 truthTables.
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    numGates = gateCounter(newCircuit)

                    if (i==1 or notAllowed.count(truthTable)==0):
                        if random.random()<0.00001:
                            print len(foundCircuits)
                        if i==1:
                            notAllowed.append(truthTable)
                        if (foundCircuits.has_key(truthTable)):
                            if numGates< gateCounter(foundCircuits[truthTable]):
                                foundCircuits[truthTable] = newCircuit
                        if (not(foundCircuits.has_key(truthTable))):
                            foundCircuits[truthTable] = newCircuit
                        if i!=5:
                            currentLevel.append(newCircuit)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                numGates = gateCounter(newCircuit)

                if (i==1 or notAllowed.count(truthTable)==0):
                    if random.random()<0.00001:
                        print len(foundCircuits)
                    if i==1:
                        notAllowed.append(truthTable)
                    if (foundCircuits.has_key(truthTable)):
                        if numGates< gateCounter(foundCircuits[truthTable]):
                            foundCircuits[truthTable] = newCircuit
                    if (not(foundCircuits.has_key(truthTable))):
                        foundCircuits[truthTable] = newCircuit
                    if i!=5:
                        currentLevel.append(newCircuit)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

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
    circuit = circuit.replace("a","00001111")
    circuit = circuit.replace("b","00110011")
    circuit = circuit.replace("c","01010101")
    while len(circuit) != 8:
        if circuit.find(truthTable)>0:
            return True
        indeces = findSmallestParentheses(circuit)
        toNor = circuit[indeces[0]:indeces[1]+1]
        answer = nor(toNor)
        circuit = circuit.replace(toNor,answer)
    return False

def levelMaker7(numLevels):
    #does not allow things with the same truthValues as levels 1 and 2.
    #if something has the same truthValues as one of its components, don't include it.
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
##    for i in range(256):
##        foundCircuits[toBinary(i)] = []
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    numGates = gateCounter(newCircuit)

                    if ((i==1 or notAllowed.count(truthTable)==0) and (not containsEqualSubCircuit(newCircuit))):
                        if random.random()<0.00001:
                            print str(len(currentLevel)) + " elements in this level so far"
                            print str(len(foundCircuits)) + " truth sets found so far"
                            print "Combining element " + str(indexAlpha+1) + " with item " +str(indexBeta+1)+ " in level " + str(j+1)
                        if i==1:
                            notAllowed.append(truthTable)
                        if (foundCircuits.has_key(truthTable)):
                            if numGates< gateCounter(foundCircuits[truthTable]):
                                foundCircuits[truthTable] = newCircuit
                        if (not(foundCircuits.has_key(truthTable))):
                            foundCircuits[truthTable] = newCircuit
                        currentLevel.append(newCircuit)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                numGates = gateCounter(newCircuit)

                if ((i==1 or notAllowed.count(truthTable)==0) and (not containsEqualSubCircuit(newCircuit))):
                    if random.random()<0.00001:
                        print str(len(currentLevel)) + " elements in this level so far"
                        print str(len(foundCircuits)) + " truth sets found so far"
                        print "Combining element " + str(indexAlpha+1) + " with item " +str(indexBeta+1)+ " in level " + str(i)
                    if i==1:
                        notAllowed.append(truthTable)
                    if (foundCircuits.has_key(truthTable)):
                        if numGates< gateCounter(foundCircuits[truthTable]):
                            foundCircuits[truthTable] = newCircuit
                    if (not(foundCircuits.has_key(truthTable))):
                        foundCircuits[truthTable] = newCircuit
                    currentLevel.append(newCircuit)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker8(numLevels, numOperations):
##    this will not save level 6 but will evaluate it up to a certain number
##    of operations and show how much time it will take
##    this will not use a dictionary. without narrowing down as we go
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    evaluations = 0
    startTime = time.time()
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    if i!=5:
                        currentLevel.append("("+alphaCircuit+"."+betaCircuit+")")
                    if i==5:
                        evaluations += 1
                        if evaluations>=numOperations:
                            endTime = time.time()
                            totalTime = endTime-startTime
                            return str(totalTime) + " seconds"
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                if i!=5:
                    currentLevel.append("("+alphaCircuit+"."+betaCircuit+")")
                if i==5:
                    evaluations += 1
                    if evaluations>=numOperations:
                        endTime = time.time()
                        totalTime = endTime-startTime
                        return str(totalTime) + " seconds"
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return allLevels

def levelMaker9(numLevels, numOperations):
##    this will not save level 6 but will evaluate it up to a certain number
##    of operations and show how much time it will take
##    this will use a dictionary. without narrowing down as we go
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    evaluations = 0
    startTime = time.time()
    for i in xrange(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in xrange(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in xrange(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in xrange(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
##                    truthTable = getTruthValue(newCircuit)
##                    if (not foundCircuits.has_key(truthTable)):
##                        foundCircuits[truthTable] = newCircuit
##                    if foundCircuits.has_key(truthTable):
##                        if gateCounter(foundCircuits[truthTable])>gateCounter(newCircuit):
##                            foundCircuits[truthTable] = newCircuit
                    if i!=5:
                        currentLevel.append(newCircuit)
                    if i==5:
                        evaluations += 1
                        if evaluations>=numOperations:
                            endTime = time.time()
                            totalTime = endTime-startTime
                            return str(totalTime) + " seconds"
            levelLookingAt = allLevels[i-1]
            for indexBeta in xrange(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
##                if (not foundCircuits.has_key(truthTable)):
##                    foundCircuits[truthTable] = newCircuit
##                if foundCircuits.has_key(truthTable):
##                    if gateCounter(foundCircuits[truthTable])>gateCounter(newCircuit):
##                        foundCircuits[truthTable] = newCircuit
                if i!=5:
                    currentLevel.append(newCircuit)
                if i==5:
                    evaluations += 1
                    if evaluations>=numOperations:
                        endTime = time.time()
                        totalTime = endTime-startTime
                        return str(totalTime) + " seconds"
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker10(numLevels, numOperations):
##    this will not save level 6 but will evaluate it up to a certain number
##    of operations and show how much time it will take
##    this will not use a dictionary. with narrowing down as we go method 1 and 4
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    evaluations = 0
    startTime = time.time()
    notAllowed = [getTruthValue('0'),getTruthValue('1'),
                  getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue('(a.b)'),getTruthValue('(a.c)'),getTruthValue('(a.0)'),
                  getTruthValue('(b.c)'),getTruthValue('(b.0)'),getTruthValue('(c.0)')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    if ((i==1)or (notAllowed.count(truthTable)==0 and not(containsEqualSubCircuit(newCircuit)))):
                        if i!=5:
                            currentLevel.append(newCircuit)
                        if i==5:
                            evaluations += 1
                            if evaluations>=numOperations:
                                endTime = time.time()
                                totalTime = endTime-startTime
                                return str(totalTime) + " seconds"
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                if ((i==1)or (notAllowed.count(truthTable)==0 and not(containsEqualSubCircuit(newCircuit)))):
                    if i!=5:
                        currentLevel.append(newCircuit)
                    if i==5:
                        evaluations += 1
                        if evaluations>=numOperations:
                            endTime = time.time()
                            totalTime = endTime-startTime
                            return str(totalTime) + " seconds"
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return allLevels

def levelMaker11(numLevels, numOperations):
##    this will not save level 6 but will evaluate it up to a certain number
##    of operations and show how much time it will take
##    this will use a dictionary. with narrowing down as we go method 1 and 4
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('0'),getTruthValue('1'),
                  getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue('(a.b)'),getTruthValue('(a.c)'),getTruthValue('(a.0)'),
                  getTruthValue('(b.c)'),getTruthValue('(b.0)'),getTruthValue('(c.0)')]
    evaluations = 0
    startTime = time.time()
    for i in xrange(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in xrange(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in xrange(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in xrange(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    if ((i==1)or (notAllowed.count(truthTable)==0 and not(containsEqualSubCircuit(newCircuit)))):
                        if (not foundCircuits.has_key(truthTable)):
                            foundCircuits[truthTable] = newCircuit
                        if foundCircuits.has_key(truthTable):
                            if gateCounter(foundCircuits[truthTable])>gateCounter(newCircuit):
                                foundCircuits[truthTable] = newCircuit
                        if i!=5:
                            currentLevel.append(newCircuit)
                        if i==5:
                            evaluations += 1
                            if evaluations>=numOperations:
                                endTime = time.time()
                                totalTime = endTime-startTime
                                return str(totalTime) + " seconds"
            levelLookingAt = allLevels[i-1]
            for indexBeta in xrange(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                if ((i==1)or (notAllowed.count(truthTable)==0 and not(containsEqualSubCircuit(newCircuit)))):
                    if (not foundCircuits.has_key(truthTable)):
                        foundCircuits[truthTable] = newCircuit
                    if foundCircuits.has_key(truthTable):
                        if gateCounter(foundCircuits[truthTable])>gateCounter(newCircuit):
                            foundCircuits[truthTable] = newCircuit
                    if i!=5:
                        currentLevel.append(newCircuit)
                    if i==5:
                        evaluations += 1
                        if evaluations>=numOperations:
                            endTime = time.time()
                            totalTime = endTime-startTime
                            return str(totalTime) + " seconds"
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def levelMaker12(num):
    allLevels = {}
    currentLevel = ['a','b','c','0']
    allLevels[0] = currentLevel
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('0'),getTruthValue('1'),
                  getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue('(a.b)'),getTruthValue('(a.c)'),getTruthValue('(a.0)'),
                  getTruthValue('(b.c)'),getTruthValue('(b.0)'),getTruthValue('(c.0)')]
    startTime = time.time()
    levelGates=0
    while not(levelGates>=num+1 or len(foundCircuits)>=256):
        
        currentLevel = allLevels[levelGates]
        print len(currentLevel)
        for alphaCircuit in currentLevel:
            alphaTruthTable = getTruthValue(alphaCircuit)
            if not(alphaTruthTable in foundCircuits):
                foundCircuits[alphaTruthTable] = alphaCircuit
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
                    for betaCircuit in comparingLevel:
                        if betaCircuit == alphaCircuit:
                            break
                        newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                        numGates = gateCounter(newCircuit)
                        if (numGates>num or (numGates>1 and(containsEqualSubCircuit(newCircuit) or (getTruthValue(newCircuit) in notAllowed)))):
                            continue
                        
                        try:
                            allLevels[numGates].append(newCircuit)
                        except KeyError:
                            allLevels[numGates] = [newCircuit]
                    break
                elif comparingLevel!=currentLevel:
                    for betaCircuit in comparingLevel:
                        newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                        numGates = gateCounter(newCircuit)
                        if (numGates>num or (numGates>1 and(containsEqualSubCircuit(newCircuit) or (getTruthValue(newCircuit) in notAllowed)))):
                            continue
                        
                        try:
                            allLevels[numGates].append(newCircuit)
                        except KeyError:
                            allLevels[numGates] = [newCircuit]
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

            
def sortAndPrint(numGates):
    foundCircuits = levelMaker(numGates)
    truths = foundCircuits.keys()
    truths.sort()
    for truth in truths:
        print truth,foundCircuits[truth]
def sortAndPrint2(foundCircuits):
    truths = foundCircuits.keys()
    truths.sort()
    for truth in truths:
        print truth,foundCircuits[truth]
    
def levelMaker13(numLevels):
    #levelMaker2 but only adds it to the list of circuits if it was useful in making the dictionary
    allLevels = []
    currentLevel = ['a','b','c','0']
    allLevels.append(currentLevel)
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),getTruthValue('0'),getTruthValue('1')]
    for i in range(1,numLevels):
        currentLevel = []
        lastLevel = allLevels[i-1]
        for indexAlpha in range(len(lastLevel)):
            alphaCircuit = lastLevel[indexAlpha]
            for j in range(0,i-1):
                levelLookingAt = allLevels[j]
                for indexBeta in range(len(levelLookingAt)):
                    betaCircuit = levelLookingAt[indexBeta]
                    newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                    truthTable = getTruthValue(newCircuit)
                    numGates = gateCounter(newCircuit)

                    if (i==1 or notAllowed.count(truthTable)==0):
                        if i==1:
                            notAllowed.append(truthTable)
                        if (foundCircuits.has_key(truthTable)):
                            if numGates< gateCounter(foundCircuits[truthTable]):
                                foundCircuits[truthTable] = newCircuit
                                currentLevel.append(newCircuit)
                        if (not(foundCircuits.has_key(truthTable))):
                            foundCircuits[truthTable] = newCircuit
                            currentLevel.append(newCircuit)
            levelLookingAt = allLevels[i-1]
            for indexBeta in range(indexAlpha+1,len(lastLevel)):
                betaCircuit = levelLookingAt[indexBeta]
                newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                truthTable = getTruthValue(newCircuit)
                numGates = gateCounter(newCircuit)

                if (i==1 or notAllowed.count(truthTable)==0):
                    if i==1:
                        notAllowed.append(truthTable)
                    if (foundCircuits.has_key(truthTable)):
                        if numGates< gateCounter(foundCircuits[truthTable]):
                            foundCircuits[truthTable] = newCircuit
                            currentLevel.append(newCircuit)
                    if (not(foundCircuits.has_key(truthTable))):
                        foundCircuits[truthTable] = newCircuit
                        currentLevel.append(newCircuit)
        allLevels.append(currentLevel)
        print "Done with level "+ str(i+1)
        print str(len(currentLevel)) + " elements at this level"
    return foundCircuits

def get228FromFile(fileName):
    foundCircuits = {}
    myFile = open(fileName,'r')
    for line in myFile:
        foundCircuits[line[0:8]] = line[9:-1]
    print len(foundCircuits)
    return foundCircuits

def checkForUse(gateCount,notUsed):
    allFoundCircuits = levelMaker15(gateCount,notUsed)[1]
    truthCircuits = get228FromFile('C:/Users/Arinze/SkyDrive/UROP_Summer_2015/GateMinimization/Stock/7gates.txt')
    unused = []
    used = []
    printed = []
    print "start"
    for circuit in allFoundCircuits[gateCount]:
        for usedCircuit in truthCircuits.values():
            if usedCircuit.count(circuit)>0:
                if used.count(circuit)==0:
##                    print "part of 228"
##                    print usedCircuit
##                    print "used"
##                    print circuit
##                    print "other similar circuits"
##                    for circuit2 in allFoundCircuits[gateCount]:
##                        if getTruthValue(circuit2)==getTruthValue(circuit) and circuit2 != circuit:
##                            print circuit2
##                            printed.append(circuit2)
##                    print "stop"
                    used.append(circuit)
                    break
    for circuit in allFoundCircuits[gateCount]:
        if used.count(circuit)==0:
            unused.append(circuit)
##    for circuit in allFoundCircuits[gateCount]:
##        if printed.count(circuit)==0 and used.count(circuit)==0:
##            print circuit
    print "Of the "+str(len(allFoundCircuits[gateCount]))+" total cicuits with "+str(gateCount)+" gates,"
    print str(len(unused))+" unused circuits with "+str(gateCount)+" gates"
    return unused#,used,printed

def levelMaker15(num, notUsed):
    allLevels = {}
    currentLevel = ['a','b','c','0']
    allLevels[0] = currentLevel
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = 'a'
    foundCircuits[getTruthValue('b')] = 'b'
    foundCircuits[getTruthValue('c')] = 'c'
    foundCircuits[getTruthValue('0')] = '0'
    foundCircuits[getTruthValue('1')] = '1'
    notAllowed = [getTruthValue('0'),getTruthValue('1'),
                  getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue('(a.b)'),getTruthValue('(a.c)'),getTruthValue('(a.0)'),
                  getTruthValue('(b.c)'),getTruthValue('(b.0)'),getTruthValue('(c.0)')]
    startTime = time.time()
    levelGates=0
    while not(levelGates>=num+1 or len(foundCircuits)>=256):
        
        currentLevel = allLevels[levelGates]
        for alphaCircuit in currentLevel:
            alphaTruthTable = getTruthValue(alphaCircuit)
            if not(alphaTruthTable in foundCircuits):
                foundCircuits[alphaTruthTable] = alphaCircuit
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
                    for betaCircuit in comparingLevel:
                        if betaCircuit == alphaCircuit:
                            break
                        newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                        numGates = gateCounter(newCircuit)
                        if ((newCircuit in notUsed) or numGates>num or (numGates>1 and(containsEqualSubCircuit(newCircuit) or (getTruthValue(newCircuit) in notAllowed)))):
                            continue
                        
                        try:
                            allLevels[numGates].append(newCircuit)
                        except KeyError:
                            allLevels[numGates] = [newCircuit]
                    break
                elif comparingLevel!=currentLevel:
                    for betaCircuit in comparingLevel:
                        newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                        numGates = gateCounter(newCircuit)
                        if ((newCircuit in notUsed) or numGates>num or (numGates>1 and(containsEqualSubCircuit(newCircuit) or (getTruthValue(newCircuit) in notAllowed)))):
                            continue
                        
                        try:
                            allLevels[numGates].append(newCircuit)
                        except KeyError:
                            allLevels[numGates] = [newCircuit]
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

def checkForUse2(gateCount,notUsed):
    allFoundCircuits = levelMaker15(gateCount,notUsed)[1]
    truthCircuits = get228FromFile('C:/Users/Arinze/SkyDrive/UROP_Summer_2015/GateMinimization/WhatToKeep/7gates.txt')
    unused = []
    used = []
    printed = []
    print "start"
    for circuit in allFoundCircuits[gateCount]:
        for usedCircuit in truthCircuits.values():
            if usedCircuit.count(circuit)>0:
                if used.count(circuit)==0:
##                    print "part of 228"
##                    print usedCircuit
##                    print "used"
##                    print circuit
##                    print "other similar circuits"
##                    for circuit2 in allFoundCircuits[gateCount]:
##                        if getTruthValue(circuit2)==getTruthValue(circuit) and circuit2 != circuit:
##                            print circuit2
##                            printed.append(circuit2)
##                    print "stop"
                    used.append(circuit)
                    break
    for circuit in allFoundCircuits[gateCount]:
        if used.count(circuit)==0:
            unused.append(circuit)
##    for circuit in allFoundCircuits[gateCount]:
##        if printed.count(circuit)==0 and used.count(circuit)==0:
##            print circuit
    print "Of the "+str(len(allFoundCircuits[gateCount]))+" total cicuits with "+str(gateCount)+" gates,"
    print str(len(unused))+" unused circuits with "+str(gateCount)+" gates"
    return unused,used#,printed
