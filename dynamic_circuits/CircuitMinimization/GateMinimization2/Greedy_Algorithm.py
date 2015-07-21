import random
import time
import cProfile
import pstats
import copy
import itertools

#Helper functions for counting gates and checking truth value
def getTruthValue(circuit):
    circuit = circuit.replace("0","00000000")
    circuit = circuit.replace("1","11111111")
    circuit = circuit.replace("a","00001111")
    circuit = circuit.replace("b","00110011")
    circuit = circuit.replace("c","01010101")
    while len(circuit) != 8:
##        print circuit 
        indeces = findSmallestParentheses(circuit)
        toNor = circuit[indeces[0]:indeces[1]+1]
        answer = nor(toNor)
        circuit = circuit.replace(toNor,answer)
    return circuit

def gateCounter(circuit):
##    return circuit.count(".")
    indeces = findLeftParentheses(circuit)
    tokens = []
    while indeces != None:
        if tokens.count(circuit[indeces[0]:indeces[1]+1]) == 0:
            tokens.append(circuit[indeces[0]:indeces[1]+1])
        circuit = circuit[0:indeces[0]] + circuit[indeces[0]+1:indeces[1]]+circuit[indeces[1]+1:len(circuit)]
        
        indeces = findLeftParentheses(circuit)
    return len(tokens)

def gateCounter2(circuit):
##    return circuit.count(".")
    indeces = findSmallestParentheses(circuit)
    count = 1
    if indeces == None:
        return 0
    while indeces != None:
##        print circuit
        circuit = circuit.replace(circuit[indeces[0]:indeces[1]+1],str(count))
        indeces = findSmallestParentheses(circuit)
        count += 1
##    print circuit
    return int(circuit)

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
##def findSmallestParentheses2(circuit):
##    #Assuming valid pairing of parentheses
##    L = len(circuit)
##    endIndex = circuit.find(")")
##    if endIndex == -1:
##        return None
##    for i in xrange(L):
##        if circuit[endIndex-i]=="(":
##            return (endIndex-i,endIndex)
##    return None
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

def sortAndPrint(foundCircuits):
    truths = foundCircuits.keys()
    truths.sort()
    for truth in truths:
        print truth," : ",foundCircuits[truth]
        
###This may need modification
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
        quoteIndeces = findLeftSingleQuotes(unformatedValue)
        while quoteIndeces != None:
            #we dont want to include the quotation marks themselves
            value.append(unformatedValue[quoteIndeces[0]+1:quoteIndeces[1]])
            unformatedValue = unformatedValue[quoteIndeces[1]+1:len(unformatedValue)]
            quoteIndeces = findLeftSingleQuotes(unformatedValue)
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

#minCircuitFinders
def minCircuitFinder5(desiredTruthValue, num=5):
    allLevels = {}
    currentLevel = ['a','b','c','0']
    allLevels[0] = currentLevel
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = ['a']
    foundCircuits[getTruthValue('b')] = ['b']
    foundCircuits[getTruthValue('c')] = ['c']
    foundCircuits[getTruthValue('0')] = ['0']
    foundCircuits[getTruthValue('1')] = ['1']
    notAllowed = [getTruthValue('0'),getTruthValue('1'),
                  getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue('(a.b)'),getTruthValue('(a.c)'),getTruthValue('(a.0)'),
                  getTruthValue('(b.c)'),getTruthValue('(b.0)'),getTruthValue('(c.0)')]
##    startTime = time.time()
    levelGates=0
    while not(levelGates>=num+1 or len(foundCircuits)>=256):
        foundTruths = []
        currentLevel = allLevels[levelGates]
        for alphaCircuit in currentLevel:
            alphaTruthValue = getTruthValue(alphaCircuit)
            if not(alphaTruthValue in foundCircuits):
                foundCircuits[alphaTruthValue] = [alphaCircuit]
                foundTruths.append(alphaTruthValue)
            elif alphaTruthValue in foundCircuits:
                if alphaTruthValue in foundTruths:
                    foundCircuits[alphaTruthValue].append(alphaCircuit)
        if len(foundCircuits)>=256:
            return foundCircuits#,allLevels
        if desiredTruthValue in foundCircuits:
            print "found it."
            return foundCircuits[desiredTruthValue]
        if levelGates == num:
##            endTime = time.time()
##            totalTime = endTime-startTime
##            print "It took "+str(totalTime)+" to finish."
            print "couldn't find it."
            return foundCircuits#,allLevels
        print "Found " + str(len(foundCircuits))+" truthValues found so far."
        for alphaCircuit in currentLevel:
            for comparingLevelGates in range(levelGates+1):
                comparingLevel = allLevels[comparingLevelGates]
                if comparingLevel==currentLevel:
                    for betaCircuitIndex in range(comparingLevel.index(alphaCircuit)+1,len(comparingLevel)):
                        betaCircuit = comparingLevel[betaCircuitIndex]
                        if betaCircuit == alphaCircuit:
                            break
                        newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                        numGates = gateCounter2(newCircuit)
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
                        numGates = gateCounter2(newCircuit)
                        if (numGates>num or (numGates>1 and(containsEqualSubCircuit(newCircuit) or (getTruthValue(newCircuit) in notAllowed)))):
                            continue
                        
                        try:
                            allLevels[numGates].append(newCircuit)
                        except KeyError:
                            allLevels[numGates] = [newCircuit]
##        print "Done with "+str(levelGates)+" gate circuits"
##        print str(levelGates+1)+" has "+str(len(allLevels[levelGates+1]))+" circuits."
##        endTime = time.time()
##        totalTime = endTime-startTime
##        print "It took "+str(totalTime)+" to get this far."

        levelGates += 1
    return foundCircuits#, allLevels
##main = []
##tried = []
def workBackward(truthValue, availableCircuits = ["00000000","11111111","00001111","00110011","01010101"], notAllowed = [],circuit = []):
    #we do not want to follow a path that leads us in a cycle. Look out for repeats.
    if truthValue in availableCircuits:
##        print "done"
        return [truthValue]
##    if truthValue in notAllowed:
##        return [None]
##    if truthValue in main:
##        return [None]
    notAllowedHere = copy.deepcopy(notAllowed)
    notAllowedHere.append(truthValue)
##    main.append(truthValue)


##    print truthValue
##    print "not allowed here"
##    print notAllowedHere
    print len(notAllowedHere)
    pairs = possiblePairs(truthValue, notAllowedHere)
    circuit = copy.deepcopy(pairs)
##    print circuit
##    pause = raw_input("pause")
    circuit2 = []
    for pair in circuit:
##        print sub
        for truth in pair:
##            print subsub
##            if subsub in tried:
##                continue
##            tried.append(subsub)
            if truth in notAllowedHere:
                break
            if type(truth)==str:
                answer = workBackward(truth,availableCircuits=availableCircuits,notAllowed=notAllowedHere)
                for pair2 in circuit:
                    for truth2 in pair2:
                        if truth2 == truth:
                            pair[pair2.index(truth)] = copy.deepcopy(answer)
            
            
##            if answer[0]==None:
##                if subsub not in notAllowedHere:
##                    notAllowedHere.append(subsub)
##                break
##        if len(sub[0])==1 and len(sub[1])==1:
##            print "Almost"
##            if sub[0][0]==None or sub[1][0]==None:
##                print "None"
##    print "End"
##    if len(circuit2)==0:
##        return [None]
    return circuit
def replaceAtPoint(s, index, new):
    return s[:index]+new+s[index+1:]
def hasX(pairs):
    for pair in pairs:
        if pair[0].count("x")>0:
            return True,pairs.index(pair)
    return [False,-1]
def possiblePairs(truthValue, notAllowedHere):
    #will print out the possible pairs that can give rise to that nor value.
    base = ""
    for char in truthValue:
        #can only get 1 if both terms are 0
        if char == "1":
            base += "0"
        if char == "0":
            base += "x"
    x = base.count("x")
    pairs = []
    pairs.append([base,base])
    while hasX(pairs)[0]:
##        print pairs
##        print hasX(pairs)
        ind = hasX(pairs)[1]
        pair = pairs[ind][:]
##        print pair[0]
        tempBase = pair[0]
        index = tempBase.find("x")
##        print tempBase
##        print index
        composite00 = replaceAtPoint(pair[0], index,"0")
        composite01 = replaceAtPoint(pair[0], index,"1")
        composite10 = replaceAtPoint(pair[1], index,"0")
        composite11 = replaceAtPoint(pair[1], index,"1")
        #at least 1 has to have a 1 at that position
        temp1 = [composite00, composite11]
        temp2 = [composite01, composite10]
        temp3 = [composite01, composite11]
        pairs.append(temp1)
        pairs.append(temp2)
        pairs.append(temp3)
        
        pairs.remove(pair)
    pairs2 = copy.deepcopy(pairs)
    for pair in pairs:
##        print pair
        if (pair[0] in notAllowedHere) or (pair[1] in notAllowedHere):
##            print "here"
            pairs2.remove(pair)
    pairs = copy.deepcopy(pairs2)
    for pair in pairs:
        reverse = [pair[1],pair[0]]
##        print pair
        if (reverse in pairs2) or (pair[0] == pair[1]):
##            print "here"
            try:
                pairs2.remove(pair)
##                print "here1"
            except:
                pass
        
##        print pairs
##        pause = raw_input("pause")
##    print len(pairs)
##    print pairs
##    if len(pairs)>0:
##        pairs = [copy.deepcopy(pairs[0])]
    return pairs2

def workBackward2(truthValue, availableCircuits = ["0000","1111","0011","0101"], notAllowed = [],circuit = []):
    #we do not want to follow a path that leads us in a cycle. Look out for repeats.
    if truthValue in availableCircuits:
        return truthValue
    notAllowedHere = copy.deepcopy(notAllowed)
    notAllowedHere.append(truthValue)
    
##    print len(notAllowedHere)
    
    pairs = possiblePairs(truthValue, notAllowedHere)
    if len(pairs)==0:
        return truthValue
    circuit = copy.deepcopy(pairs)
    circuit2 = copy.deepcopy(circuit)
    for pair in circuit:
        for truth in pair:
            if truth in notAllowedHere:
                break
            if type(truth)==str:
                answer = workBackward2(truth,availableCircuits=availableCircuits,notAllowed=notAllowedHere)
                if type(answer) != str:
                    for pair2 in circuit:
                        for truth2 in pair2:
                            if truth2 == truth:
                                pair[pair2.index(truth)] = copy.deepcopy(answer)
                elif type(answer)==str:
                    if answer not in availableCircuits:
                        try:
                            circuit2.remove(pair)
                            break
                        except:
                            break
    circuit = copy.deepcopy(circuit2)
                        
    if len(circuit)==0:
        return truthValue
    return circuit

def workBackward3(truthValue, availableCircuits = ["0000","1111","0011","0101"], notAllowed = []):
    #made it work now I need to figure out how to parse the list of lists returned and find the smallest circuit.
    if truthValue in availableCircuits:
##        print "done"
        return truthValue
    if truthValue in notAllowed:
##        print "backin' up"
        return truthValue
    notAllowed = notAllowed + [truthValue]
    replacement = possiblePairs(truthValue,notAllowed)
    if len(replacement)==0:
##        print "backin' up 2"
        return truthValue
    for indexPair in range(len(replacement)):
        pair = replacement[indexPair]
        for indexTruth in range(2):
            truth = pair[indexTruth]
            answer = workBackward3(truth,availableCircuits=availableCircuits,notAllowed=notAllowed)
            pair[indexTruth] = answer
    replacement2 = copy.deepcopy(replacement)
    for indexPair in range(len(replacement2)):
        pair = replacement2[indexPair]
        for indexTruth in range(2):
            truth = pair[indexTruth]
            if type(truth) == str:
                if truth not in availableCircuits:
                    replacement.remove(pair)
                    #may give an error
                    break

    if len(replacement)==0:
##        print "backin' up 3"
        return truthValue
    for indexPair in range(len(replacement)):
        replacement[indexPair] = tuple(replacement[indexPair])
    return replacement
            
def workBackward4(truthValue,depthLimit, availableCircuits = ["00000000","11111111","00001111","00110011","01010101"], notAllowed = [],depth=0):
    #workBackward3 for 3 inputs
    depth += 1
    if depth >= depthLimit:
        return truthValue
    if truthValue in availableCircuits:
        return truthValue
    if truthValue in notAllowed:
        return truthValue
    notAllowed = notAllowed + [truthValue]
    replacement = possiblePairs(truthValue,notAllowed)
    if len(replacement)==0:
        return truthValue
    for indexPair in range(len(replacement)):
        pair = replacement[indexPair]
        for indexTruth in range(2):
            truth = pair[indexTruth]
            answer = workBackward4(truth,depthLimit,availableCircuits=availableCircuits,notAllowed=notAllowed, depth=depth)
            pair[indexTruth] = answer
    replacement2 = copy.deepcopy(replacement)
    for indexPair in range(len(replacement2)):
        pair = replacement2[indexPair]
        for indexTruth in range(2):
            truth = pair[indexTruth]
            if type(truth) == str:
                if truth not in availableCircuits:
                    replacement.remove(pair)
                    #may give an error
                    break

    if len(replacement)==0:
        return truthValue
    for indexPair in range(len(replacement)):
        replacement[indexPair] = tuple(replacement[indexPair])
    return replacement
def combine(truthValue, depthLimit=3,num=5,askToStop=False):
    """
        suggested depthLimit is 3.
        suggested num is 5.
    """
    startTime = time.time()
    x = minCircuitFinder5(truthValue, num=num)
    if type(x) == list:
        print "done"
        return x
    elif type(x) == dict:
        x = minCircuitFinder6(truthValue, num=num)
    print "starting advanced algorithm"
    availableCircuits = x.keys()
    for i in range(3,depthLimit+1):
        unformatted = workBackward4(truthValue,i, availableCircuits=availableCircuits)
        if unformatted != truthValue:
            break
    if unformatted == truthValue:
        print "the depthLimit was not great enough"
        return
    print "checkpoint 1"
    formatted = reformat2(unformatted)
    minCircuits = []
    minNumGates = 999
    oldMinNumGates = minNumGates
    print "checkpoint 2"
    for i in xrange(len(formatted)):
        
        if i%500==0 or (i<500 and i%50==0):
            print str(i),"out of",str(len(formatted))
            print "so far the best circuit has",str(minNumGates),"gates"
            oldMinNumGates = minNumGates
            if askToStop:
                currentTime = time.time()
                print "This has been running for",str(currentTime - startTime),"seconds"
                if minNumGates == num+1:
                    answer = ""
                    while (answer != "Y" and answer != "N"):
                        answer = raw_input("We have found at least one circuit with the minimum number of gates.\nDo you want to continue and find all the circuits with this number of gates?\n Y=Yes\nN=No\nYour answer: ")
                    if answer == "N":
                        return minCircuits
                    elif answer=="Y":
                        oldMinNumGates = minNumGates

                elif((currentTime - startTime >= 180) and (oldMinNumGates == minNumGates)) :
                    #if it has been running for more than 3 minutes ask to stop
                    answer = ""
                    while (answer != "1" and answer != "2"):
                        answer = raw_input("This has run for more than 3 minutes and the gate count has not reduced. This may be the lowest number of gate we can find.\nDo you want to continue and find the minimum still, or stop now?\n1=continue\n2=stop now\nYour answer: ")
                    if answer == "2":
                        return minCircuits
                    elif answer=="1":
                        oldMinNumGates = minNumGates
                    
        s = formatted[i]
        circuits,gateCount = makeMinCircuit(s,x)
        if gateCount<minNumGates:
            minCircuits = circuits
            minNumGates = gateCount
        elif gateCount==minNumGates:
            minCircuits = minCircuits + circuits
    if len(minCircuits)==0:
        print "something isn't right"
    print minNumGates
    print "done2"
    endTime = time.time()
    print str(endTime - startTime),"seconds"
    return minCircuits



##def convertToCircuit(s,circuitDict,minNum=10):
##    startS = ""
##    copyS = s
##    while (startS != s):
##        startS = s
##        length = len(s)
##        index = length-1
##        
##        while index>=0:
##            if s[index]==")":
##                if s[index-23]=="(":
##                    temp = s[index-23:index]
##                    circ1 = circuitDict[temp[2:10]][0]
##                    circ2 = circuitDict[temp[14:22]][0]
##                    built = "("+circ1+"."+circ2+")"
##    ##                print built
##                    numGates = gateCounter(built)
##                    if numGates>=minNum:
##    ##                    print s[index-23:index+4]
##    ##                    pause = raw_input("pause")
##                        try:
##                            s = s[:index-23]+s[index+3:]
##                        except:
##                            s = s[:index-23]+s[index+1:]
##    ##                    print gateCounter(built)
##                        index = index-23
##                    elif numGates<minNum:
##                        minNum = numGates
####                        s = s[:index-23]+built+s[index+1:]
##                        index = index-23
##            index -=1
##    s = copyS
####    print minNum
##    index = len(s)-1
##    while index>=0:
##        if s[index]==")":
##            if s[index-23]=="(":
##                temp = s[index-23:index]
##                circ1 = circuitDict[temp[2:10]][0]
##                circ2 = circuitDict[temp[14:22]][0]
##                built = "("+circ1+"."+circ2+")"
##                numGates = gateCounter(built)
####                print numGates
##                if numGates>minNum:
##                    if s[index+1]==",":
##                        s = s[:index-23]+s[index+3:]
##                    else:
##                        s = s[:index-23]+s[index+1:]
##                    index = index-23
##                elif numGates<=minNum:
##                    minNum -=1
##                    s = s[:index-23]+built+s[index+1:]
####                    return built ,minNum  
##        index -=1
##    return s, minNum
            
def findLeftBracket(s):
    L = len(s)
    counter = 0
    indeces = []
    for i in xrange(L):
        if s[i] == "[":
            counter += 1
            if len(indeces)==0:
                indeces.append(i)
        elif s[i] == "]":
            counter -= 1
            if counter == 0:
                indeces.append(i)
                return indeces
    return None    
        
def reformat(data):
    copyData = str(data)
    results = [copyData]
    results2 = copy.deepcopy(results)
##    x=0
    resultsStart = []
    while resultsStart != results:
        resultsStart = copy.deepcopy(results)
##        print results
##        pause = raw_input("pause1")
        for item in results:
##            print item
##            pause = raw_input("pause2")
            copyItem = item
            bracketIndeces = findLeftBracket(copyItem)
            if bracketIndeces != None:
##                print "here"
                choices = []
                toAdd = []
                inBrackets = copyItem[bracketIndeces[0]:bracketIndeces[1]+1]
                inBrackets2 = inBrackets
                parenIndeces = findLeftParentheses(inBrackets)
                while parenIndeces != None:
                    choices.append(inBrackets[parenIndeces[0]:parenIndeces[1]+1])
                    inBrackets = inBrackets[0:parenIndeces[0]] + inBrackets[parenIndeces[1]+1:len(inBrackets)]
                    parenIndeces = findLeftParentheses(inBrackets)
                for choice in choices:
                    toAdd.append(copyItem.replace(inBrackets2,choice))
                results2.remove(item)
                for newItem in toAdd:
                    results2.append(newItem)
##                print results2
##                pause = raw_input("pause3")
        results = copy.deepcopy(results2)
    return results

def reformat2(data):
    copyData = str(data)
    results = [copyData]
##    results2 = copy.deepcopy(results)
##    x=0
##    resultsStart = []
    i=0
    while i<len(results):
        item = results[i]
        i+=1
        copyItem = item
        bracketIndeces = findLeftBracket(copyItem)
        if bracketIndeces != None:
            i-=1
            choices = []
            toAdd = []
            inBrackets = copyItem[bracketIndeces[0]:bracketIndeces[1]+1]
            inBrackets2 = inBrackets
            parenIndeces = findLeftParentheses(inBrackets)
            while parenIndeces != None:
                choices.append(inBrackets[parenIndeces[0]:parenIndeces[1]+1])
                inBrackets = inBrackets[0:parenIndeces[0]] + inBrackets[parenIndeces[1]+1:len(inBrackets)]
                parenIndeces = findLeftParentheses(inBrackets)
            for choice in choices:
                toAdd.append(copyItem.replace(inBrackets2,choice))
            results.remove(item)
            for newItem in toAdd:
                results.append(newItem)
    return results

def findSmallestApostrophe(s):
    #Assuming valid pairing of parentheses
    L = len(s)
    startIndex = -1
    for i in xrange(L):
        if s[i]=="'":
            startIndex = i
            break
    if startIndex == -1:
        return None
    for i in xrange(startIndex+1,L):
        if s[i]=="'":
            endIndex = i
            return (startIndex,endIndex)
    return None

def makeMinCircuit(s,foundCircuits):
    s1 = s
    truths = []
    potentialCircuits = []
    apIndeces = findSmallestApostrophe(s1)
    minCircuits = []
    minNumGates = 999
    while apIndeces != None:
        tempTruth = s1[apIndeces[0]+1:apIndeces[1]]
        if tempTruth not in truths:
            truths.append(tempTruth)
        s1 = s1[0:apIndeces[0]] + s1[apIndeces[1]+1:len(s1)]
        apIndeces = findSmallestApostrophe(s1)
##    print truths
    for truth in truths:
        potentialCircuits.append(foundCircuits[truth])
    combinations = list(itertools.product(*potentialCircuits))
    s = s.replace(", ",".")
    s = s.replace("'","")
##    print s
    for combo in combinations:
##        pause = raw_input("pause1")
        s3 = s
        for i in range(len(combo)):
            s3 = s3.replace(truths[i],combo[i])
##        print s3
        numGates = gateCounter(s3)
##        print numGates
        if numGates < minNumGates:
            minCircuits = [s3]
            minNumGates = numGates
        elif numGates == minNumGates:
            minCircuits.append(s3)
    if len(minCircuits)==0:
        print "something happened"
    return minCircuits, minNumGates

def combine2(truthValue, depthLimit=3,num=5,askToStop=False, foundCircuits = {}):
    """
        suggested depthLimit is 3.
        suggested num is 5.
    """
    startTime = time.time()
    if foundCircuits == {}:
        x = minCircuitFinder5(truthValue, num=num)
    elif foundCircuits!={}:
        x = foundCircuits
    if type(x) == list:
##        print "done"
        return x
##    print "starting advanced algorithm"
    availableCircuits = x.keys()
    for i in range(depthLimit+1):
        unformatted = workBackward4(truthValue,i, availableCircuits=availableCircuits)
        if unformatted != truthValue:
            break
    if unformatted == truthValue:
##        print "the depthLimit was not great enough"
        return
##    print "checkpoint 1"
    formatted = reformat2(unformatted)
    minCircuits = []
    minNumGates = 999
    oldMinNumGates = minNumGates
##    print "checkpoint 2"
    count = 0
    for i in xrange(len(formatted)):
##        if i%10==0:
##            print str(i),"out of",str(len(formatted))
##            print "so far the best circuit has",str(minNumGates),"gates"
##            if (oldMinNumGates == minNumGates):
##                count+=1
##            elif (oldMinNumGates != minNumGates):
##                oldMinNumGates = minNumGates
##                count = 0
        if askToStop:
##                currentTime = time.time()
##                print i, currentTime - startTime
##                print "This has been running for",str(currentTime - startTime),"seconds"
##                if count>=7:
##                    return minCircuits
            if minNumGates == num+1:
                answer = "N"
                while (answer != "Y" and answer != "N"):
                    answer = raw_input("We have found at least one circuit with the minimum number of gates.\nDo you want to continue and find all the circuits with this number of gates?\n Y=Yes\nN=No\nYour answer: ")
                if answer == "N":
                    return minCircuits
                elif answer=="Y":
                    oldMinNumGates = minNumGates

            elif (time.time() - startTime >= 200):
                #if it has been running for more than 200 seconds ask to stop
                answer = "2"
                while (answer != "1" and answer != "2"):
                    answer = raw_input("This has run for more than 3 minutes and the gate count has not reduced. This may be the lowest number of gate we can find.\nDo you want to continue and find the minimum still, or stop now?\n1=continue\n2=stop now\nYour answer: ")
                if answer == "2":
                    return minCircuits
                elif answer=="1":
                    oldMinNumGates = minNumGates
                    
        s = formatted[i]
        circuits,gateCount = makeMinCircuit(s,x)
        if gateCount<minNumGates:
            minCircuits = circuits
            minNumGates = gateCount
        elif gateCount==minNumGates:
            minCircuits = minCircuits + circuits
##    if len(minCircuits)==0:
##        print "something isn't right"
##    print minNumGates
##    print "done2"
##    endTime = time.time()
##    print str(endTime - startTime),"seconds"
    return minCircuits        

def minCircuitFinder6(desiredTruthValue, num=5):
    allLevels = {}
    currentLevel = ['a','b','c','0']
    allLevels[0] = currentLevel
    foundCircuits = {}
    foundCircuits[getTruthValue('a')] = ['a']
    foundCircuits[getTruthValue('b')] = ['b']
    foundCircuits[getTruthValue('c')] = ['c']
    foundCircuits[getTruthValue('0')] = ['0']
    foundCircuits[getTruthValue('1')] = ['1']
    notAllowed = [getTruthValue('0'),getTruthValue('1'),
                  getTruthValue('a'),getTruthValue('b'),getTruthValue('c'),
                  getTruthValue('(a.b)'),getTruthValue('(a.c)'),getTruthValue('(a.0)'),
                  getTruthValue('(b.c)'),getTruthValue('(b.0)'),getTruthValue('(c.0)')]
##    startTime = time.time()
    levelGates=0
    while not(levelGates>=num+1 or len(foundCircuits)>=256):
##        foundTruths = []
        currentLevel = allLevels[levelGates]
        for alphaCircuit in currentLevel:
            alphaTruthValue = getTruthValue(alphaCircuit)
            if not(alphaTruthValue in foundCircuits):
                foundCircuits[alphaTruthValue] = [alphaCircuit]
##                foundTruths.append(alphaTruthValue)
            elif alphaTruthValue in foundCircuits:
##                if alphaTruthValue in foundTruths:
                foundCircuits[alphaTruthValue].append(alphaCircuit)
        if len(foundCircuits)>=256:
            return foundCircuits#,allLevels
        if desiredTruthValue in foundCircuits:
            print "found it."
            return foundCircuits[desiredTruthValue]
        if levelGates == num:
##            endTime = time.time()
##            totalTime = endTime-startTime
##            print "It took "+str(totalTime)+" to finish."
            print "couldn't find it."
            return foundCircuits#,allLevels
        print "Found " + str(len(foundCircuits))+" truthValues found so far."
        for alphaCircuit in currentLevel:
            for comparingLevelGates in range(levelGates+1):
                comparingLevel = allLevels[comparingLevelGates]
                if comparingLevel==currentLevel:
                    for betaCircuitIndex in range(comparingLevel.index(alphaCircuit)+1,len(comparingLevel)):
                        betaCircuit = comparingLevel[betaCircuitIndex]
                        if betaCircuit == alphaCircuit:
                            break
                        newCircuit = "("+alphaCircuit+"."+betaCircuit+")"
                        numGates = gateCounter2(newCircuit)
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
                        numGates = gateCounter2(newCircuit)
                        if (numGates>num or (numGates>1 and(containsEqualSubCircuit(newCircuit) or (getTruthValue(newCircuit) in notAllowed)))):
                            continue
                        
                        try:
                            allLevels[numGates].append(newCircuit)
                        except KeyError:
                            allLevels[numGates] = [newCircuit]
##        print "Done with "+str(levelGates)+" gate circuits"
##        print str(levelGates+1)+" has "+str(len(allLevels[levelGates+1]))+" circuits."
##        endTime = time.time()
##        totalTime = endTime-startTime
##        print "It took "+str(totalTime)+" to get this far."

        levelGates += 1
    return foundCircuits#, allLevels
def combine3(truthValue, depthLimit=3,num=5,askToStop=False, foundCircuits = None):
    #uses minCircuitFinder6 instead of 5
    """
        suggested depthLimit is 3.
        suggested num is 5.
    """
    startTime = time.time()
    if foundCircuits == None:
        x = minCircuitFinder6(truthValue, num=num)
    elif foundCircuits != None:
        x = foundCircuits
    if type(x) == list:
        print "done"
        return x
    print "starting advanced algorithm"
    availableCircuits = x.keys()
    for i in range(depthLimit+1):
        unformatted = workBackward4(truthValue,i, availableCircuits=availableCircuits)
        if unformatted != truthValue:
            break
    if unformatted == truthValue:
        print "the depthLimit was not great enough"
        return
    print "checkpoint 1"
    formatted = reformat2(unformatted)
    minCircuits = []
    minNumGates = 999
    oldMinNumGates = minNumGates
    print "checkpoint 2"
    for i in xrange(len(formatted)):
        
        if i%500==0 or (i<500 and i%50==0):
            print str(i),"out of",str(len(formatted))
            print "so far the best circuit has",str(minNumGates),"gates"
            oldMinNumGates = minNumGates
            if askToStop:
                currentTime = time.time()
                print "This has been running for",str(currentTime - startTime),"seconds"
                if minNumGates == num+1:
                    answer = ""
                    while (answer != "Y" and answer != "N"):
                        answer = raw_input("We have found at least one circuit with the minimum number of gates.\nDo you want to continue and find all the circuits with this number of gates?\n Y=Yes\nN=No\nYour answer: ")
                    if answer == "N":
                        return minCircuits
                    elif answer=="Y":
                        oldMinNumGates = minNumGates

                elif((currentTime - startTime >= 180) and (oldMinNumGates == minNumGates)) :
                    #if it has been running for more than 3 minutes ask to stop
                    answer = ""
                    while (answer != "1" and answer != "2"):
                        answer = raw_input("This has run for more than 3 minutes and the gate count has not reduced. This may be the lowest number of gate we can find.\nDo you want to continue and find the minimum still, or stop now?\n1=continue\n2=stop now\nYour answer: ")
                    if answer == "2":
                        return minCircuits
                    elif answer=="1":
                        oldMinNumGates = minNumGates
                    
        s = formatted[i]
        circuits,gateCount = makeMinCircuit(s,x)
        if gateCount<minNumGates:
            minCircuits = circuits
            minNumGates = gateCount
        elif gateCount==minNumGates:
            minCircuits = minCircuits + circuits
    if len(minCircuits)==0:
        print "something isn't right"
    print minNumGates
    print "done2"
    endTime = time.time()
    print str(endTime - startTime),"seconds"
    return minCircuits

def combine4(truthValue, depthLimit=3,num=5,askToStop=False):
    """
        suggested depthLimit is 3.
        suggested num is 5.
    """
    startTime = time.time()
    x = minCircuitFinder5(truthValue, num=num)
    if type(x) == list:
        print "done"
        return x
    elif type(x) == dict:
        x = minCircuitFinder6(truthValue, num=num)
    print "starting advanced algorithm"
    availableCircuits = x.keys()
    for i in range(depthLimit+1):
        unformatted = workBackward4(truthValue,i, availableCircuits=availableCircuits)
        if unformatted != truthValue:
            break
    if unformatted == truthValue:
        print "the depthLimit was not great enough"
        return
    print "checkpoint 1"
    formatted = reformat2(unformatted)
    minCircuits = []
    minNumGates = 999
    oldMinNumGates = minNumGates
    print "checkpoint 2"
    for i in xrange(len(formatted)):
        
        if i%500==0 or (i<500 and i%50==0):
            print str(i),"out of",str(len(formatted))
            print "so far the best circuit has",str(minNumGates),"gates"
            oldMinNumGates = minNumGates
            if askToStop:
                currentTime = time.time()
                print "This has been running for",str(currentTime - startTime),"seconds"
                if minNumGates == num+1:
                    answer = ""
                    while (answer != "Y" and answer != "N"):
                        answer = raw_input("We have found at least one circuit with the minimum number of gates.\nDo you want to continue and find all the circuits with this number of gates?\n Y=Yes\nN=No\nYour answer: ")
                    if answer == "N":
                        return minCircuits
                    elif answer=="Y":
                        oldMinNumGates = minNumGates

                elif((currentTime - startTime >= 180) and (oldMinNumGates == minNumGates)) :
                    #if it has been running for more than 3 minutes ask to stop
                    answer = ""
                    while (answer != "1" and answer != "2"):
                        answer = raw_input("This has run for more than 3 minutes and the gate count has not reduced. This may be the lowest number of gate we can find.\nDo you want to continue and find the minimum still, or stop now?\n1=continue\n2=stop now\nYour answer: ")
                    if answer == "2":
                        return minCircuits
                    elif answer=="1":
                        oldMinNumGates = minNumGates
                    
        s = formatted[i]
        circuits,gateCount = makeMinCircuit(s,x)
        if gateCount<minNumGates:
            minCircuits = circuits
            minNumGates = gateCount
        elif gateCount==minNumGates:
            minCircuits = minCircuits + circuits
    if len(minCircuits)==0:
        print "something isn't right"
    print minNumGates
    print "done2"
    endTime = time.time()
    print str(endTime - startTime),"seconds"
    return minCircuits
