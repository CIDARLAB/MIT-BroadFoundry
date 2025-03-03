# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 15:54:35 2015

@author: Arinze
Edited by Alex (7.29.2015)
"""

"""
This file was used to do a number of things. Its primary function is to compare
the netlists given to me by Bryan with the circuits we generated. The circuits 
are simplified by replacing the final NOT with a OR where applicable. This also 
contains a function that can convert the circuitString into a netlist, but it 
only works for NOR and OR. Minor modification required to recognize other gate 
types.
"""

import json
import Truths_and_Gates as tg

fileName = "JsonFiles/Circuits_MaxCost9.0_NorOnly.json"
placeToSaveOR = "JsonFiles/outputs/withOR.json"
placeToSaveNOR = "JsonFiles/outputs/withNOR.json"

def getCircuitsFromFile(fileName):
    myFile = open(fileName,'r')
    data = json.load(myFile)
    data = data[2:255]
    # Maps truth values to circuits
    circuits = {}

    # Add the 0 and 1 circuits
    circuits['00000000'] = ['0']
    circuits['11111111'] = ['1']
    
    # Maps 253 of the remaining 254 truth values between '0' and '1' to corresponding circuit
    for d in data:
        circuits[d["truthValue"]] = d["circuits"]

    # Maps last truth value and append additional circuits to truth values. 11101001 is the 
    # missing truth value. (Useful for an 'or' simplification in a later step)
    for d in data:
        tv = d["truthValue"]
        newtv = invert(tv)
        tempOldCirc = circuits[tv]
        tempNewCirc = []
        for circ in tempOldCirc:
            tempNewCirc.append("("+circ+".0)")
        if newtv in circuits:
            for circ in tempNewCirc:
                if circ not in circuits[newtv] and ".0).0)" not in circ:
                    circuits[newtv].append(circ)
        else:
            circuits[newtv] = tempNewCirc
    myFile.close()
    return circuits

# Function for the inversion of a truth value string    
def invert(truthValue):
    s = ""
    for char in truthValue:
        if char=="1":
            s+="0"
        else:
            s+="1"
    return s

# Function returns simplified circuit with an 'or' in place of an ending 'nor' and 'not' 
def simplify(fileName, placeToSaveOR):
    # 'or' (+) cost is free because it is just 2 promoters before the FP output
    opcost = {"~":1,"&":1,"@":1,"+":0,"^":1,".":1,"=":1,">":1,"$":1}
    circuits = getCircuitsFromFile(fileName)
    simplifiedCircuits = {}
    truthValues = circuits.keys()
    truthValues.remove("00000000")
    truthValues.remove("11111111")
    simplifiedCircuits['00000000'] = ['0']
    simplifiedCircuits['11111111'] = ['1']
    
    #replace final 'not' (.0) with 'or' (+) when applicable
    for tv in truthValues:
        currCircuits = circuits[tv]
        tempNewCircuits = []
        for circ in currCircuits:
            tempNewCircuits.append(convertToORIfPossible(circ))
        simplifiedCircuits[tv] = tempNewCircuits
        
    for tv in truthValues:
        currCircuits = simplifiedCircuits[tv]
        tempNewCircuits = [currCircuits[0]]
        currBestCost = tg.circuitCost(currCircuits[0],cost=opcost)
        for circ in currCircuits[1:]:
            cost = tg.circuitCost(circ,cost=opcost)
            if cost<currBestCost:
                currBestCost = cost
                tempNewCircuits = [circ]
            elif cost==currBestCost:
                tempNewCircuits.append(circ)
        simplifiedCircuits[tv] = tempNewCircuits
    
    finalResults = {}
    for tv in simplifiedCircuits:
        info = {}
        info["circuits"] = simplifiedCircuits[tv]
        info["cost"] = tg.circuitCost(simplifiedCircuits[tv][0],cost=opcost)
        finalResults[tv] = info
    myFile = open(placeToSaveOR,'w')
    json.dump(finalResults, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()
    return finalResults
    
def convertToORIfPossible(circuitString):
    #ignore circuits like (a.0) but not ((a.b).0)
    if circuitString[-4:] != ").0)":
        return circuitString
    #remove trailing 'not' (.0)
    circuitString = circuitString[1:-3]
    count = 0
    for i in range(1,len(circuitString)):
        char = circuitString[i]
        if char == "(":
            count += 1
        elif char == ")":
            count -= 1
        if count == 0:
            return circuitString[:i+1] + "+" +circuitString[i+2:]
    
# Adds a cost value to a new Json file. 
def remakeFile(fileName,placeToSaveNOR):
    myFile = open(fileName,'r')
    data = json.load(myFile)
    data = data[2:255]
    #maps truth values to circuits
    circuits = {}
    #add the 0 and 1 circuits
    infoOne = {}
    infoZero = {}
    infoOne["circuits"] = ['1']
    infoOne["cost"] = 0
    infoZero["circuits"] = ['0']
    infoZero["cost"] = 0
    circuits['00000000'] = infoZero
    circuits['11111111'] = infoOne
    
    # Creates dictionary for circuits with additional cost reference
    for d in data:
        tv = d["truthValue"]
        currCircuits = d["circuits"]
        if tv == "00010110":
            newCircuits = []
            for circ in currCircuits:
                newCircuits.append("("+circ+".0)")
            tempInfo = {}
            tempInfo["circuits"] = newCircuits
            tempInfo["cost"] = tg.circuitCost(newCircuits[0])
            circuits["11101001"] = tempInfo
        tempInfo = {}
        tempInfo["circuits"] = currCircuits
        tempInfo["cost"] = tg.circuitCost(currCircuits[0])
        circuits[tv] = tempInfo
    myFile.close()
    
    myFile = open(placeToSaveNOR,'w')
    json.dump(circuits, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()  
    return circuits

myFileNameNOR = "JsonFiles/withNOR.json"
bryanFileNameNOR = "JsonFiles/netlist_in3out1.json"
myFileNameOR = "JsonFiles/withOR.json"
bryanFileNameOR = "JsonFiles/netlist_in3out1_OR.json"

def compare(myFileName, bryanFileName):
    """
    takes in two file names and compares the cost of each truth value
    """
    myFile = open(myFileName,'r')
    bryanFile = open(bryanFileName,'r')
    myData = json.load(myFile)
    bryanData = json.load(bryanFile)
    
    truthValues = myData.keys()
    truthValues.sort()
    truthValues = truthValues[1:255]
    differences = {}
    for i in range(len(truthValues)):
        tv = truthValues[i]
        myCost = myData[tv]["cost"]
        bryanCost = len(bryanData[i]["netlist"])
        if bryanData[i]["netlist"][-1][0:2]=="OR":
            bryanCost -= 1
        if myData[tv]["circuits"][0] in ["a","b","c"]:
            myCost = 1
        if myCost != bryanCost:
            print "Difference in",tv,"myCost:",myCost,"bryanCost:", bryanCost
            differences[tv] = [myCost,bryanCost]
        else:
            print "truthValue:",tv,"cost:",myCost
    myFile.close()
    bryanFile.close()
    return differences
    
def convertToNetlist(circuitString):
    """
    Takes in a circuitString and converts it to a netlist. The format of the 
    netlist will be a list of strings. Only works for 'nor', 'not', and 'or'.
    """
    wireCount = 0
    netlist = []
    parenIndex = tg.findSmallestParentheses(circuitString)
    while parenIndex != None:
        wireName = "W"+str(wireCount)
        subCircuit = circuitString[parenIndex[0]:parenIndex[1]+1]
        gate = ""
        if subCircuit.find("+")!=-1:
            gate += "OR("
            operator = "+"
        elif subCircuit.find(".")!=-1:
            operator = "."
            if subCircuit.find(".0")!=-1:
                gate += "NOT("
            else:
                gate += "NOR("
        gate = gate + wireName + ","
        #remove parentheses and split by operator
        fanIn = subCircuit[1:-1]
        fanIn = fanIn.split(operator)
        
        for w in fanIn:
            if w!="0":
                gate = gate + w + ","
        #remove trailing comma and close parentheses
        gate = gate[:-1] + ")"
        
        netlist.append(gate)
        circuitString = circuitString.replace(subCircuit, wireName)
        parenIndex = tg.findSmallestParentheses(circuitString)
        wireCount += 1
        
    #Incase we are just looking at the inputs
    if wireCount == 0:
        netlist.append("BUF(y,"+circuitString+")")
    return netlist

unstandarizedNOR = "JsonFiles/withNOR.json"
unstandarizedOR = "JsonFiles/withOR.json"
standarizedNOR = "JsonFiles/standardNOR.json"
standarizedOR = "JsonFiles/standardOR.json"

def standardizeFileWithNetlists(fileName, fileDestination):
    myFile = open(fileName,'r')
    myData = json.load(myFile)
    myFile.close()
    allTruthValues = myData.keys()
    newData = {}
    
    for tv in allTruthValues:
        allCircuits = myData[tv]["circuits"]
        cost = myData[tv]["cost"]
        tempInfo = {}
        listOfNetlist = []
        for circ in allCircuits:
            listOfNetlist.append(convertToNetlist(circ))
        tempInfo["cost"] = cost
        tempInfo["netlists"] = listOfNetlist
        newData[tv] = tempInfo
    
    myFile = open(fileDestination,'w')
    json.dump(newData, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()