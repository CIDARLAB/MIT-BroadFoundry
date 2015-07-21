from Truths_and_Gates import *
from GateMinimizationIO import *
import time
import cProfile
import pstats
import copy

##I need to figure out how to generalize for number of inputs
def minCircuitFinder1(numInputs, operators= [], cost = {}, maxCost=5):
    """
        Takes the number of inputs, which operators to use, the cost of each
        operator, and the maximum cost to check up to. Default is to use no
        operators. If no weight is given for an operator, it assumes a cost of 1.
        If no maxCost is specified, 
        Returns the dictionary of truth values mapped to their minimal cost
        circuits.
    """
    startTime = time.time()
    #Calculate how many gates we are expecting to see
    tempValue = 1
    for i in range(numInputs):
        tempValue *= 2
    expectedNumTruths = 1
    for i in range(tempValue):
        expectedNumTruths *= 2
    print "We expect to see",expectedNumTruths,"possible truths values"
    
    #Stop if one of the given operators is not recognized
    approvedOperators = ["~","&","@","+","^",".","=",">","$"]
    for operator in operators:
        if operator not in approvedOperators:
            print operator, "is not a recognized operator."
            return
    #Default cost to 1 for unspecified operators or for operators whose costs
    #are negative.
    for operator in operators:
        if operator not in cost:
            cost[operator] = 1
        elif operator in cost:
            if cost[operator]<0:
                cost[operator] = 1
    #Determines whether or not we can check the circuits at the beginning or the end
    has0Penalties = False
    for operator in operators:
        if cost[operator] == 0:
            has0Penalties=True
            break
    #Remove unnecessary operators from cost dictionary
    for operator in approvedOperators:
        if operator not in operators:
            cost.pop(operator,None)
    #If nor costs 0 we can do ((a.b).c) then keep nor ing infinitely until
    #we find everything. This may not finish because it will try to enumerate
    #everything possible before stopping.
    "Check what I can make 0"
    for operator in cost.keys():
        if cost[operator]==0 and operator!="~":
            print "You cannot assign a cost of 0 to the", operator,"operator."
            return
    #Determine which of the operators used are symmetric e.g. (a+b) == (b+a)
    symmetricOp = []
    asymmetricOp = []
    specialOp = []
    for operator in operators:
        if operator in ["&","@","+","^",".","="]:
            symmetricOp.append(operator)
        elif operator in [">","$"]:
            asymmetricOp.append(operator)
        else: #operator == "~"
            specialOp.append(operator)
            
    #Initialize the dictionary with the minimum circuits with the inputs, 0, and 1
    foundTruthValues = {}
    foundTruthValues[getTruthValue("1",numInputs)] = "1"
    foundTruthValues[getTruthValue("0",numInputs)] = "0"
    foundTruthValues[getTruthValue("a",numInputs)] = "a"
    if numInputs>=2:
        foundTruthValues[getTruthValue("b",numInputs)] = "b"
    if numInputs>=3:
        foundTruthValues[getTruthValue("c",numInputs)] = "c"
    
    sortedByCost = {}
    sortedByCost[0] = ["a"]
    if numInputs>=2:
        sortedByCost[0].append("b")
    if numInputs>=3:
        sortedByCost[0].append("c")
    #Don't save anything that is just 1 or 0 or one of the inputs.
##    reallyNotAllowed = []
    reallyNotAllowed = [getTruthValue("1",numInputs),getTruthValue("0",numInputs)]#,getTruthValue("a",numInputs),getTruthValue("b",numInputs),getTruthValue("c",numInputs)]
##    #If something equals the same truth value as the smallest non-zero costing circuits...
##    I'm not sure if this is still allowed
##    notAllowed = []

    #We are going to be iterating through a list as we mutate it.
    allCosts = sortedByCost.keys()
    allCosts.sort()
    for alphaGroupCost in allCosts:
##        print allCosts
        alphaGroup = sortedByCost[alphaGroupCost]
        indexAlphaGroupCost = allCosts.index(alphaGroupCost)
        if not(has0Penalties):
            foundThisRound = []
            if indexAlphaGroupCost <= 1:
                for alphaCircuit in alphaGroup:
                    truthValue = getTruthValue(alphaCircuit,numInputs)
                    if truthValue not in reallyNotAllowed:
                        reallyNotAllowed.append(truthValue)
            for alphaCircuit in alphaGroup:
                truthValue = getTruthValue(alphaCircuit,numInputs)
                if (truthValue in foundTruthValues) and (truthValue in foundThisRound):
                    foundTruthValues[truthValue].append(alphaCircuit)
                elif (truthValue not in foundTruthValues):
                    foundTruthValues[truthValue] = [alphaCircuit]
                    foundThisRound.append(truthValue)
            currTime = time.time()
            print "----------------------------"
            print "Cost:",alphaGroupCost,"Time:",round(currTime-startTime,3),"seconds"
            print len(foundTruthValues),"truth values found so far"
            print len(alphaGroup), "circuits in this group"
            if len(foundTruthValues) == expectedNumTruths or (alphaGroupCost==maxCost):
                "Done"
                return foundTruthValues, sortedByCost

        for alphaCircuit in alphaGroup:
            #make the NOT and special considerations
            specialNewCircuits = []
            for operator in specialOp:
                specialNewCircuits.append("("+operator+alphaCircuit+")")
            for operator in symmetricOp:
                if operator == "." or operator == "=":
                    specialNewCircuits.append("("+alphaCircuit+operator+"0)")
                elif operator == "^":
                    specialNewCircuits.append("("+alphaCircuit+operator+"1)")
            for operator in asymmetricOp:
                if operator == ">":
                    specialNewCircuits.append("("+alphaCircuit+">0)")
                elif operator == "$":
                    specialNewCircuits.append("(1$"+alphaCircuit+")")
##            print specialNewCircuits
            for newCircuit in specialNewCircuits:
                newCircuitCost = circuitCost(newCircuit,cost=cost)
                newCircuitTruthValue = getTruthValue(newCircuit, numInputs)
##                print sortedByCost
##                print newCircuit
                if newCircuitCost in sortedByCost:
                    if (newCircuit not in sortedByCost[newCircuitCost]) and (newCircuitTruthValue not in reallyNotAllowed) and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                        sortedByCost[newCircuitCost].append(newCircuit)
                elif (newCircuitCost not in sortedByCost):
                    if (newCircuitTruthValue not in reallyNotAllowed) and newCircuitCost<=maxCost and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                        sortedByCost[newCircuitCost] = [newCircuit]
                        #mutate the list we are iterating through
                        allCosts.append(newCircuitCost)
                        allCosts.sort()
##                print sortedByCost
##                print newCircuit
            for indexBetaGroupCost in range(indexAlphaGroupCost+1):
                betaGroupCost = allCosts[indexBetaGroupCost]
                betaGroup = sortedByCost[betaGroupCost]
                if indexBetaGroupCost != indexAlphaGroupCost:
                    for betaCircuit in betaGroup:
                        #make all possible new circuits that don't use NOT
                        newCircuits = []
                        for operator in symmetricOp:
                            newCircuits.append("("+alphaCircuit+operator+betaCircuit+")")
                        for operator in asymmetricOp:
                            newCircuits.append("("+alphaCircuit+operator+betaCircuit+")")
                            newCircuits.append("("+betaCircuit+operator+alphaCircuit+")")
                        #add the new circuits to their correct position
##                        print newCircuits
                        for newCircuit in newCircuits:
                            newCircuitCost = circuitCost(newCircuit,cost=cost)
                            newCircuitTruthValue = getTruthValue(newCircuit, numInputs)
                            if newCircuitCost in sortedByCost:
##                                print sortedByCost
##                                print newCircuit
                                if (newCircuit not in sortedByCost[newCircuitCost]) and (newCircuitTruthValue not in reallyNotAllowed) and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost].append(newCircuit)
                            elif (newCircuitCost not in sortedByCost):
                                if (newCircuitTruthValue not in reallyNotAllowed) and newCircuitCost<=maxCost and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost] = [newCircuit]
                                    #mutate the list we are iterating through
                                    allCosts.append(newCircuitCost)
                                    allCosts.sort()
                        
                if indexBetaGroupCost == indexAlphaGroupCost:
                    #must go combine with everything to the left so we don't miss
                    #something if an operator cost is 0
                    for betaCircuit in betaGroup:
                        if alphaCircuit == betaCircuit:
                            break
                        #make all possible new circuits that don't use NOT
                        newCircuits = []
                        for operator in symmetricOp:
                            newCircuits.append("("+betaCircuit+operator+alphaCircuit+")")
                        for operator in asymmetricOp:
                            newCircuits.append("("+alphaCircuit+operator+betaCircuit+")")
                            newCircuits.append("("+betaCircuit+operator+alphaCircuit+")")
                        #add the new circuits to their correct position
##                        print newCircuits
                        for newCircuit in newCircuits:
##                            print sortedByCost
##                            print circuit
                            newCircuitCost = circuitCost(newCircuit,cost=cost)
                            newCircuitTruthValue = getTruthValue(newCircuit, numInputs)
                            if newCircuitCost in sortedByCost:
##                                print sortedByCost
##                                print newCircuit
                                if (newCircuit not in sortedByCost[newCircuitCost]) and (newCircuitTruthValue not in reallyNotAllowed) and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost].append(newCircuit)
                            elif (newCircuitCost not in sortedByCost):
                                if (newCircuitTruthValue not in reallyNotAllowed) and newCircuitCost<=maxCost and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost] = [newCircuit]
                                    #mutate the list we are iterating through
                                    allCosts.append(newCircuitCost)
                                    allCosts.sort()
        #We must check at the end to avoid missing circuits because of 0 cost operator
        if has0Penalties:
            foundThisRound = []
            if indexAlphaGroupCost <= 1:
                for alphaCircuit in alphaGroup:
                    truthValue = getTruthValue(alphaCircuit,numInputs)
                    if truthValue not in reallyNotAllowed:
                        reallyNotAllowed.append(truthValue)
            for alphaCircuit in alphaGroup:
                truthValue = getTruthValue(alphaCircuit,numInputs)
                if (truthValue in foundTruthValues) and (truthValue in foundThisRound):
                    foundTruthValues[truthValue].append(alphaCircuit)
                elif (truthValue not in foundTruthValues):
                    foundTruthValues[truthValue] = [alphaCircuit]
                    foundThisRound.append(truthValue)
##            print alphaGroup
            currTime = time.time()
            print "----------------------------"
            print "Cost:",alphaGroupCost,"Time:",round(currTime-startTime,3),"seconds"
            print len(foundTruthValues),"truth values found so far"
            print len(alphaGroup), "circuits in this group"
##            try:
##                print len(sortedByCost[allCosts[indexAlphaGroupCost+1]]), "in the next group"
##            except:
##                pass
##            print allCosts
            if len(foundTruthValues) == expectedNumTruths or (alphaGroupCost==maxCost):
                "Done"
                return foundTruthValues,sortedByCost
        

def minCircuitFinder11(numInputs, operators= [], cost = {}, maxCost=5):
    """
        Takes the number of inputs, which operators to use, the cost of each
        operator, and the maximum cost to check up to. Default is to use no
        operators. If no weight is given for an operator, it assumes a cost of 1.
        If no maxCost is specified, 
        Returns the dictionary of truth values mapped to their minimal cost
        circuits.
    """
    startTime = time.time()
    #Calculate how many gates we are expecting to see
    expectedNumTruths = 1
    for i in range(numInputs):
        expectedNumTruths *= 2
    tempValue = expectedNumTruths
    expectedNumTruths = 1
    for i in range(tempValue):
        expectedNumTruths *= 2
    print "We expect to see",expectedNumTruths,"possible truths values"
    
    #Stop if one of the given operators is not recognized
    approvedOperators = ["~","&","@","+","^",".","=",">","$"]
    for operator in operators:
        if operator not in approvedOperators:
            print operator, "is not a recognized operator."
            return
    #Default cost to 1 for unspecified operators or for operators whose costs
    #are negative.
    for operator in operators:
        if operator not in cost:
            cost[operator] = 1
        elif operator in cost:
            if cost[operator]<0:
                cost[operator] = 1
    #Determines whether or not we can check the circuits at the beginning or the end
    has0Penalties = False
    for operator in operators:
        if cost[operator] == 0:
            has0Penalties=True
            break
    #Remove unnecessary operators from cost dictionary
    for operator in approvedOperators:
        if operator not in operators:
            cost.pop(operator,None)
    #If nor costs 0 we can do ((a.b).c) then keep nor ing infinitely until
    #we find everything. This may not finish because it will try to enumerate
    #everything possible before stopping.
    "Check what I can make 0"
    for operator in cost.keys():
        if cost[operator]==0 and operator!="~":
            print "You cannot assign a cost of 0 to the", operator,"operator."
            return
    #Determine which of the operators used are symmetric e.g. (a+b) == (b+a)
    symmetricOp = []
    asymmetricOp = []
    specialOp = []
    for operator in operators:
        if operator in ["&","@","+","^",".","="]:
            symmetricOp.append(operator)
        elif operator in [">","$"]:
            asymmetricOp.append(operator)
        else: #operator == "~"
            specialOp.append(operator)
            
    #Initialize the dictionary with the minimum circuits with the inputs, 0, and 1
    foundTruthValues = {}
    foundTruthValues[getTruthValue("1",numInputs)] = "1"
    foundTruthValues[getTruthValue("0",numInputs)] = "0"
    foundTruthValues[getTruthValue("a",numInputs)] = "a"
    if numInputs>=2:
        foundTruthValues[getTruthValue("b",numInputs)] = "b"
    if numInputs>=3:
        foundTruthValues[getTruthValue("c",numInputs)] = "c"
    
    sortedByCost = {}
    sortedByCost[0] = ["a"]
    if numInputs>=2:
        sortedByCost[0].append("b")
    if numInputs>=3:
        sortedByCost[0].append("c")
    #Don't save anything that is just 1 or 0 or one of the inputs.
##    reallyNotAllowed = []
    reallyNotAllowed = [getTruthValue("1",numInputs),getTruthValue("0",numInputs)]#,getTruthValue("a",numInputs),getTruthValue("b",numInputs),getTruthValue("c",numInputs)]
##    #If something equals the same truth value as the smallest non-zero costing circuits...
##    I'm not sure if this is still allowed
##    notAllowed = []

    #We are going to be iterating through a list as we mutate it.
    allCosts = sortedByCost.keys()
    allCosts.sort()
    for alphaGroupCost in allCosts:
##        print allCosts
        alphaGroup = sortedByCost[alphaGroupCost]
        indexAlphaGroupCost = allCosts.index(alphaGroupCost)
        if not(has0Penalties):
            foundThisRound = []
            if indexAlphaGroupCost <= 1:
                for alphaCircuit in alphaGroup:
                    truthValue = getTruthValue(alphaCircuit,numInputs)
                    if truthValue not in reallyNotAllowed:
                        reallyNotAllowed.append(truthValue)
            alphaTruth = str(alphaGroup)
            alphaTruth = alphaTruth.replace("'","")
            alphaTruth = alphaTruth.replace("[","")
            alphaTruth = alphaTruth.replace("]","")
            alphaTruth = getTruthValue(alphaTruth,numInputs)
            alphaTruth = alphaTruth.split(", ")
            for i in range(len(alphaGroup)):
                alphaCircuit = alphaGroup[i]
                truthValue=alphaTruth[i]
                if (truthValue in foundTruthValues) and (truthValue in foundThisRound):
                    foundTruthValues[truthValue].append(alphaCircuit)
                elif (truthValue not in foundTruthValues):
                    foundTruthValues[truthValue] = [alphaCircuit]
                    foundThisRound.append(truthValue)
            currTime = time.time()
            print "----------------------------"
            print "Cost:",alphaGroupCost,"Time:",round(currTime-startTime,3),"seconds"
            print len(foundTruthValues),"truth values found so far"
            print len(alphaGroup), "circuits in this group"
            if len(foundTruthValues) == expectedNumTruths or (alphaGroupCost==maxCost):
                "Done"
                return foundTruthValues, sortedByCost

        for alphaCircuit in alphaGroup:
            #make the NOT and special considerations
            specialNewCircuits = []
            for operator in specialOp:
                specialNewCircuits.append("("+operator+alphaCircuit+")")
            for operator in symmetricOp:
                if operator == "." or operator == "=":
                    specialNewCircuits.append("("+alphaCircuit+operator+"0)")
                elif operator == "^":
                    specialNewCircuits.append("("+alphaCircuit+operator+"1)")
            for operator in asymmetricOp:
                if operator == ">":
                    specialNewCircuits.append("("+alphaCircuit+">0)")
                elif operator == "$":
                    specialNewCircuits.append("(1$"+alphaCircuit+")")
##            print specialNewCircuits
            for newCircuit in specialNewCircuits:
                newCircuitCost = circuitCost(newCircuit,cost=cost)
                newCircuitTruthValue = getTruthValue(newCircuit, numInputs)
##                print sortedByCost
##                print newCircuit
                if newCircuitCost in sortedByCost:
                    if (newCircuit not in sortedByCost[newCircuitCost]) and (newCircuitTruthValue not in reallyNotAllowed) and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                        sortedByCost[newCircuitCost].append(newCircuit)
                elif (newCircuitCost not in sortedByCost):
                    if (newCircuitTruthValue not in reallyNotAllowed) and newCircuitCost<=maxCost and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                        sortedByCost[newCircuitCost] = [newCircuit]
                        #mutate the list we are iterating through
                        allCosts.append(newCircuitCost)
                        allCosts.sort()
##                print sortedByCost
##                print newCircuit
            for indexBetaGroupCost in range(indexAlphaGroupCost+1):
                betaGroupCost = allCosts[indexBetaGroupCost]
                betaGroup = sortedByCost[betaGroupCost]
                if indexBetaGroupCost != indexAlphaGroupCost:
                    for betaCircuit in betaGroup:
                        #make all possible new circuits that don't use NOT
                        newCircuits = []
                        for operator in symmetricOp:
                            newCircuits.append("("+alphaCircuit+operator+betaCircuit+")")
                        for operator in asymmetricOp:
                            newCircuits.append("("+alphaCircuit+operator+betaCircuit+")")
                            newCircuits.append("("+betaCircuit+operator+alphaCircuit+")")
                        #add the new circuits to their correct position
##                        print newCircuits
                        for newCircuit in newCircuits:
                            newCircuitCost = circuitCost(newCircuit,cost=cost)
                            newCircuitTruthValue = getTruthValue(newCircuit, numInputs)
                            if newCircuitCost in sortedByCost:
##                                print sortedByCost
##                                print newCircuit
                                if (newCircuit not in sortedByCost[newCircuitCost]) and (newCircuitTruthValue not in reallyNotAllowed) and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost].append(newCircuit)
                            elif (newCircuitCost not in sortedByCost):
                                if (newCircuitTruthValue not in reallyNotAllowed) and newCircuitCost<=maxCost and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost] = [newCircuit]
                                    #mutate the list we are iterating through
                                    allCosts.append(newCircuitCost)
                                    allCosts.sort()
                        
                if indexBetaGroupCost == indexAlphaGroupCost:
                    #must go combine with everything to the left so we don't miss
                    #something if an operator cost is 0
                    for betaCircuit in betaGroup:
                        if alphaCircuit == betaCircuit:
                            break
                        #make all possible new circuits that don't use NOT
                        newCircuits = []
                        for operator in symmetricOp:
                            newCircuits.append("("+betaCircuit+operator+alphaCircuit+")")
                        for operator in asymmetricOp:
                            newCircuits.append("("+alphaCircuit+operator+betaCircuit+")")
                            newCircuits.append("("+betaCircuit+operator+alphaCircuit+")")
                        #add the new circuits to their correct position
##                        print newCircuits
                        for newCircuit in newCircuits:
##                            print sortedByCost
##                            print circuit
                            newCircuitCost = circuitCost(newCircuit,cost=cost)
                            newCircuitTruthValue = getTruthValue(newCircuit, numInputs)
                            if newCircuitCost in sortedByCost:
##                                print sortedByCost
##                                print newCircuit
                                if (newCircuit not in sortedByCost[newCircuitCost]) and (newCircuitTruthValue not in reallyNotAllowed) and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost].append(newCircuit)
                            elif (newCircuitCost not in sortedByCost):
                                if (newCircuitTruthValue not in reallyNotAllowed) and newCircuitCost<=maxCost and (not(containsEqualSubCircuit(newCircuit,numInputs))):
                                    sortedByCost[newCircuitCost] = [newCircuit]
                                    #mutate the list we are iterating through
                                    allCosts.append(newCircuitCost)
                                    allCosts.sort()
        #We must check at the end to avoid missing circuits because of 0 cost operator
        if has0Penalties:
            foundThisRound = []
            if indexAlphaGroupCost <= 1:
                for alphaCircuit in alphaGroup:
                    truthValue = getTruthValue(alphaCircuit,numInputs)
                    if truthValue not in reallyNotAllowed:
                        reallyNotAllowed.append(truthValue)
            for alphaCircuit in alphaGroup:
                truthValue = getTruthValue(alphaCircuit,numInputs)
                if (truthValue in foundTruthValues) and (truthValue in foundThisRound):
                    foundTruthValues[truthValue].append(alphaCircuit)
                elif (truthValue not in foundTruthValues):
                    foundTruthValues[truthValue] = [alphaCircuit]
                    foundThisRound.append(truthValue)
##            print alphaGroup
            currTime = time.time()
            print "----------------------------"
            print "Cost:",alphaGroupCost,"Time:",round(currTime-startTime,3),"seconds"
            print len(foundTruthValues),"truth values found so far"
            print len(alphaGroup), "circuits in this group"
##            try:
##                print len(sortedByCost[allCosts[indexAlphaGroupCost+1]]), "in the next group"
##            except:
##                pass
##            print allCosts
            if len(foundTruthValues) == expectedNumTruths or (alphaGroupCost==maxCost):
                "Done"
                return foundTruthValues,sortedByCost


def minCircuitFinder2():
    """
        Takes a truth value, the number of inputs, which operators to use,
        the cost of each operator, and the maximum cost to check up to.
        Default is to use no operators.
        Returns the list of circuits with the minimal cost for the truth
        value if the truth value was found. If the truth value was not found,
        returns the dictionary of truth values mapped to their minimal cost
        circuits.
    """
    pass
