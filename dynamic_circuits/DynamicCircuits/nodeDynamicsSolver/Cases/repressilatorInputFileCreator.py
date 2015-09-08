# -*- coding: utf-8 -*-
"""
Created on Friday Sept 4, 2015

NOTE: Uses nodeDynamicSolver.py to display repressilator circuit dynamics.

@author: Alex Lim
File that is used to create master input file for plotting repressilator dynamics.

Parameters came from the paper: 
A synthetic oscillatory network of transcriptional regulators (Elowitz and 
Leibler, 2000).
"""

import json
import numpy as np
import random

def writeToJson(inputObj,fileLoc):
    """
    Writes the object, preferably composed of lists or dictionaries, to a json file with pretty print
    """
    myFile = open(fileLoc,'w')
    json.dump(inputObj, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()

#Initialization of variables for genes
nodes = ['X','Y','Z']
iniCond = {}
inputD = {}

for node in nodes:
	iniCond[node] = random.uniform(0.0,300.0)			#initial condition of gene
	inputD[node] = {}
	inputD[node]['name']= node
	inputD[node]['hl'] = random.uniform(0.75,1.25)*10.0 	# half-life (minutes)
	inputD[node]['a'] = np.log(2)/inputD[node]['hl']		# degradation rate
	inputD[node]['B0'] = random.uniform(0.75,1.25)*0.6		# basal production rate
	inputD[node]['B'] = random.uniform(0.75,1.25)*300.0 	# induced production rate
	inputD[node]['n'] = random.uniform(0.75,1.25)*2.0 		# hill coefficient
	inputD[node]['Km'] = random.uniform(0.75,1.25)*40.0 	# coefficient for activation or repression

inputObject = [iniCond, inputD]

addOn = []

# Creates an individual case of the necessary node interactions. 
# In this case it is a repressilator (X-|Y-|Z-|X)

current = {}
current['from'] = nodes[0]
current['to'] = nodes[1]
current['effect'] = -1
addOn.append(current)

current = {}
current['from'] = nodes[1]
current['to'] = nodes[2]
current['effect'] = -1
addOn.append(current)

current = {}
current['from'] = nodes[2]
current['to'] = nodes[0]
current['effect'] = -1
addOn.append(current)
inputObject.append(addOn)

#Creates the master input file
fileName = 'repressilator.json'
writeToJson(inputObject,fileName)