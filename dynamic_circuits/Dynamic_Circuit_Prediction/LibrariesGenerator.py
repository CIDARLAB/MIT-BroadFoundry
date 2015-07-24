import json
import random

standardNonInput = {"Km":1000,"n":2,"mB":30,"pB":20,"m_halfLife":2,"p_halfLife":10}
standardInput = {"Km":5,"n":2,"mB":30}
fileLoc = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015"
numInputs = 4
numRepressors = 10
numOutputs = 3

def generateLibraries(numInputs,numRepressors, numOutputs, fileLoc, version, variance=5.0, standardNonInput=standardNonInput, standardInput=standardInput):
    """
    generates libraries using the standard by choosing random values within the
    specified variance. Values will be in the range of
    standard/variance<value<variance*standard
    """
    
    inputBaseName = "IPTG"
    repressorBaseName = "GENE"
    outputBaseName = "FP"
    
    inputProperties = ["Km","n","mB"]
    repressorProperties = ["Km","n","mB","pB","m_halfLife","p_halfLife"]
    outputProperties = ["pB","m_halfLife","p_halfLife"]
    
    
    inputList = []
    repressorList = []
    outputList = []
    
    varRange = [1/float(variance), float(variance)]
    
    for i in range(numInputs):
        tempInput = {}
        tempInput["NAME"] = inputBaseName+str(i+1)
        tempInput["TYPE"] = "INPUT"
        for prop in inputProperties:
            if prop != "n":
                tempInput[prop] = standardInput[prop]*random.uniform(varRange[0],varRange[1])
            else:
                tempInput[prop] = random.uniform(1,4)
        inputList.append(tempInput)
        
            
    for i in range(numRepressors):
        tempRepressor = {}
        tempRepressor["NAME"] = repressorBaseName+str(i+1)
        tempRepressor["TYPE"] = "REPRESSOR"
        for prop in repressorProperties:
            if prop != "n":
                tempRepressor[prop] = standardNonInput[prop]*random.uniform(varRange[0],varRange[1])
            else:
                tempRepressor[prop] = random.uniform(1,4)            
        repressorList.append(tempRepressor)
        
        
    for i in range(numOutputs):
        tempOutput = {}
        tempOutput["NAME"] = outputBaseName+str(i+1)
        tempOutput["TYPE"] = "OUTPUT"
        for prop in outputProperties:
            tempOutput[prop] = standardNonInput[prop]*random.uniform(varRange[0],varRange[1])
        outputList.append(tempOutput)
    
    if fileLoc[-1]!="/":
        fileLoc += "/"
    inputsFileName = fileLoc + "InputLibrary" + str(version) + ".json"
    repressorFileName = fileLoc + "RepressorLibrary" + str(version) + ".json"
    OutputFileName = fileLoc + "OutputLibrary" + str(version) + ".json"
    
    myFile1 = open(inputsFileName,'w')
    json.dump(inputList, myFile1, sort_keys=True, indent=4, separators=(',', ': '))
    myFile1.close()
    myFile2 = open(repressorFileName,'w')
    json.dump(repressorList, myFile2, sort_keys=True, indent=4, separators=(',', ': '))
    myFile2.close()
    myFile3 = open(OutputFileName,'w')
    json.dump(outputList, myFile3, sort_keys=True, indent=4, separators=(',', ': '))
    myFile3.close()