# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 12:40:24 2015

@author: Arinze
"""

"""
This file was used to take in the netlists given to me by bryan and split them
up into individual files named by their truth value and containing only their
netlist.
"""

import json


def getAllBinaries(numInputs=3):
    """
    returns a list of all binary values that can represent a truth
    value given the number of inputs excluding all 0 and all 1
    """
    binaryLength = 2**numInputs
    maxVal = 2**binaryLength
    binaries = []
    for i in range(1,maxVal-1):
        temp="{0:b}".format(i)
        while len(temp)<binaryLength:
            temp = "0" + temp
        binaries.append(temp)
    return binaries
        
        
def splitAndSave(fileDir, outputFileDir):
    """
    Takes in the netlist tha Bryan sent me and splits them into individual files
    """
    if not(outputFileDir.endswith("/") or outputFileDir.endswith("\\")):
        outputFileDir += "/"
    myFile = open(fileDir,'r')
    allData = json.load(myFile)
#    print allData[1]
    
    binaries = getAllBinaries(3)

    for i in range(len(binaries)):
        tempTruthValue = binaries[i]
        tempNetlist = allData[i]["netlist"]
        outputFile = open(outputFileDir + tempTruthValue + ".json",'w')
        json.dump(tempNetlist, outputFile, sort_keys=True, indent=4, separators=(',', ': '))

outputFileDir_OR = "JsonFiles/SplitByTruthValue_OR"
outputFileDir = "JsonFiles/SplitByTruthValue"
fileDir = "JsonFiles/netlist_in3out1.json"
fileDir_OR = "JsonFiles/netlist_in3out1_OR.json"