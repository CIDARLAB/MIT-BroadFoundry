import os
from GateMinimizationAND import *

#base text: os.system('C:/Users/Arinze/Desktop/abc10216.exe -c "read C:/Users/Arinze/Desktop/test.blif; strash; write C:/Users/Arinze/Desktop/out.bench; quit;"')

def toBinary(num):
    #converts a decimal number to an 8 bit binary string
    s = "{0:b}".format(num)
    while len(s)<8:
        s = "0"+s
    return s

## write blif file
def makeBlifFor(truthValue, directoryBlif="C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ABC/Blif"):
    """
        This will take in a truth value as a string and create a blif file for it.
        The directory cannot have folders whose names start with numbers.
    """
    abc = ["000","001","010","011","100","101","110","111"]
    directoryBlif.replace("\\","/")
    if directoryBlif[-1] != "/":
        directoryBlif += '/'
    myFile = open(directoryBlif+truthValue+".blif", "w")
    myFile.write(".model blif_file\n.inputs in1 in2 in3\n.outputs out\n.names in1 in2 in3 out\n")
    for i in range(8):
        if truthValue[i]=="1":
            myFile.write(abc[i]+" 1\n")
    myFile.write(".end")
    myFile.close()
    
def convertBlifToBench(truthValue, directoryBench="C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ABC/Bench", directoryBlif="C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ABC/Blif"):
    """
        Finds the blif file for that truth value and converts it to a bench
    """
    ABC = "C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ABC/abc10216.exe"
    directoryBench.replace("\\","/")
    directoryBlif.replace("\\","/")
    if directoryBench[-1] != "/":
        directoryBench += '/'
    if directoryBlif[-1] != "/":
        directoryBlif += '/'
    saveBenchIn = directoryBench + truthValue+".bench"
    getBlifFrom = directoryBlif + truthValue+".blif"
    os.system(ABC + ' -c "read '+getBlifFrom+'; strash; rewrite; refactor; balance; write '+saveBenchIn+'; quit;"')
    
def makeBenchForAll():
    for i in range(1,255):
        truthValue = toBinary(i)
        makeBlifFor(truthValue)
        convertBlifToBench(truthValue)

def compareGateCount():
    andGateCountsABC = {}
    andGateCountsFound = {}
    different = []
    notFound = []
    for i in range(1,255):
        truthValue = toBinary(i)
        myFile = open("C:/Users/Arinze/SkyDrive/UROP_Summer_2015/ABC/Bench/"+truthValue+".bench","r")
        numAnds = 0
        for line in myFile:
            if line.count("AND")>0:
                numAnds += 1
        andGateCountsABC[truthValue] = numAnds
    #when stock 5gates is made change file name to 5gates and also make 6 from 5 and try again with that.
    myFile = open("C:/Users/Arinze/SkyDrive/UROP_Summer_2015/GateMinimizationAND/from5to7fixed/6gates.txt","r")
    for line in myFile:
        key = line[0:8]
        rest = line[8:]
        indeces = findLeftDoubleQuotes(rest)
        if indeces == None:
            indeces = findLeftSingleQuotes(rest)
        circuit = rest[indeces[0]+1:indeces[1]]
        numAnds = gateCounter(circuit)
        andGateCountsFound[key] = numAnds
    keys = andGateCountsABC.keys()
    keys.sort()
    print "truthVal, ABC : Found"
    for key in keys:
        try:
##            print key+", "+str(andGateCountsABC[key])+" : "+str(andGateCountsFound[key])
            if andGateCountsABC[key] != andGateCountsFound[key]:
                different.append(key)
        except:
##            print key+", "+str(andGateCountsABC[key])+" : x"
            notFound.append(key)
    for key in different:
        print key+", "+str(andGateCountsABC[key])+" : "+str(andGateCountsFound[key])
    for key in notFound:
        print key+", "+str(andGateCountsABC[key])+" : x"
            
    return andGateCountsABC, andGateCountsFound
        
