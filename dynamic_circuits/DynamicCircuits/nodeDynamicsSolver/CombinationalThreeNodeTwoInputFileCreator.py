# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
File that is used to create master input file needed for XY input plots, but can be used 
as basis for converting objects to json files.

Parameters came from the paper: 
A synthetic oscillatory network of transcriptional regulators (Elowitz and 
Leibler, 2000).
"""

import json
import numpy as np
import random

# Derived from Arinze from Online
def writeToJson(inputObj,fileLoc):
    """
    Writes the object, preferably composed of lists or dictionaries, to a json file with pretty print
    """
    myFile = open(fileLoc,'w')
    json.dump(inputObj, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()

#Initialization of variables for genes, addition of boolean input value for node
nodes = ['X','Y','Z']
inputs = ['A','B']
iniCond = {}
nodeParams = {}
inputParams = {}

#Creates random parameters for the genes
for node in nodes:
	iniCond[node] = random.uniform(0.0,300.0)			#initial condition of gene
	nodeParams[node] = {}
	nodeParams[node]['name']= node
	nodeParams[node]['hl'] = random.uniform(0.75,1.25)*10.0 	# half-life (minutes)
	nodeParams[node]['a'] = np.log(2)/nodeParams[node]['hl']		# degradation rate
	nodeParams[node]['B0'] = random.uniform(0.75,1.25)*0.0		# basal production rate
	nodeParams[node]['B'] = random.uniform(0.75,1.25)*300.0 	# induced production rate
	nodeParams[node]['n'] = random.uniform(0.75,1.25)*2.0 		# hill coefficient
	nodeParams[node]['Km'] = random.uniform(0.75,1.25)*40.0 	# coefficient for activation or repression

for inp in inputs:
    inputParams[inp] = {}
    inputParams[inp]['name'] = inp
    inputParams[inp]['B'] = random.uniform(0.75,1.25)*300.0 	# induced production rate
    inputParams[inp]['n'] = random.uniform(0.75,1.25)*2.0 		# hill coefficient
    inputParams[inp]['Km'] = random.uniform(0.75,1.25)*40.0 	# coefficient for activation
    # For the signal, the value at each time point would be within a few magnitudes of Km
    if(inp == 'A'):
          inputParams[inp]['signal'] = [0,10,0,10]
    else:
          inputParams[inp]['signal'] = [0,0,10,10]
inputObject = [iniCond, nodeParams, inputParams]


# Creates an individual case of all node interactions. 
# In this case, a combinational circuit 	Xi-->A-->C
#								Yi-->B--7
inputToNode = []
nodeToNode = []

'''
#Random node vector assignment
for inputFrom in inputs:
	for nodeTo in nodes:
		inputToNode = {}
		inputToNode['from'] = inputFrom
		inputToNode['to'] = nodeTo
		current['effect'] = random.randint(0,1)  # On or Off
		nodeToNode.append(current)
'''

#input to node interactions
AtoX = 1
AtoY = 0
BtoX = 0
BtoY = 1

current = {}
current['from'] = 'A'
current['to'] = 'X'
current['effect'] = AtoX
inputToNode.append(current)

current = {}
current['from'] = 'A'
current['to'] = 'Y'
current['effect'] = AtoY
inputToNode.append(current)

current = {}
current['from'] = 'B'
current['to'] = 'X'
current['effect'] = BtoX
inputToNode.append(current)

current = {}
current['from'] = 'B'
current['to'] = 'Y'
current['effect'] = BtoY
inputToNode.append(current)

#node to node interactions
XtoX = 0
XtoY = 0
XtoZ = -1
YtoX = 0
YtoY = 0
YtoZ = -1
ZtoX = 0
ZtoY = 0
ZtoZ = 0

current = {}
current['from'] = 'X'
current['to'] = 'X'
current['effect'] = XtoX
nodeToNode.append(current)

current = {}
current['from'] = 'X'
current['to'] = 'Y'
current['effect'] = XtoY
nodeToNode.append(current)

current = {}
current['from'] = 'X'
current['to'] = 'Z'
current['effect'] = XtoZ
nodeToNode.append(current)

current = {}
current['from'] = 'Y'
current['to'] = 'X'
current['effect'] = YtoX
nodeToNode.append(current)

current = {}
current['from'] = 'Y'
current['to'] = 'Y'
current['effect'] = YtoY
nodeToNode.append(current)

current = {}
current['from'] = 'Y'
current['to'] = 'Z'
current['effect'] = YtoZ
nodeToNode.append(current)

current = {}
current['from'] = 'Z'
current['to'] = 'X'
current['effect'] = ZtoX
nodeToNode.append(current)

current = {}
current['from'] = 'Z'
current['to'] = 'Y'
current['effect'] = ZtoY
nodeToNode.append(current)

current = {}
current['from'] = 'Z'
current['to'] = 'Z'
current['effect'] = ZtoZ
nodeToNode.append(current)

inputObject.append(inputToNode)
inputObject.append(nodeToNode)

#Creates the master input file consisting of initial conditions, node then input parameters, 
# and then input to node then node to node vectors
fileName = 'masterCombiNodeInputFile.json'
writeToJson(inputObject,fileName)