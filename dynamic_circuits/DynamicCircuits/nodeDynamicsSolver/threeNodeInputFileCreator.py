# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
File that is used to create master input file needed for XY input plots, but can be used 
as basis for converting objects to json files.
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
'''
for node in nodes:
	iniCond[node] = random.uniform(0.0,10.0)				#initial condition of gene
	inputD[node] = {}
	inputD[node]['name']= node
	inputD[node]['hl'] = 10.0 	# half-life (minutes)
	inputD[node]['a'] = np.log(2)/inputD[node]['hl']		# degradation rate
	inputD[node]['B0'] = 0.6		# basal production rate
	inputD[node]['B'] = 300.0 	# induced production rate
	inputD[node]['n'] = 2.0 		# hill coefficient
	inputD[node]['Km'] = 40.0 	# coefficient for activation or repression
'''
inputObject = [iniCond, inputD]

'''
Think about metrics/characteristics such as dampening coefficient, amplitude, period
for X,Y, and Z.
Script that runs through n metrics, and then convert to jpeg, gif
'''

addOn = []
'''
#Random node vector assignment
for nodeFrom in var:
	for nodeTo in var:
		current = {}
		current['from'] = nodeFrom
		current['to'] = nodeTo
		current['effect'] = random.randint(-1,1)
		addOn.append(current) 
'''
# Creates an individual case of all node interactions. 
# In this case it is a repressilator (X-|Y-|Z-|X)
XtoX = 0
XtoY = -1
XtoZ = 0
YtoX = 0
YtoY = 0
YtoZ = -1
ZtoX = -1
ZtoY = 0
ZtoZ = 0

current = {}
current['from'] = 'X'
current['to'] = 'X'
current['effect'] = XtoX
addOn.append(current)

current = {}
current['from'] = 'X'
current['to'] = 'Y'
current['effect'] = XtoY
addOn.append(current)

current = {}
current['from'] = 'X'
current['to'] = 'Z'
current['effect'] = XtoZ
addOn.append(current)

current = {}
current['from'] = 'Y'
current['to'] = 'X'
current['effect'] = YtoX
addOn.append(current)

current = {}
current['from'] = 'Y'
current['to'] = 'Y'
current['effect'] = YtoY
addOn.append(current)

current = {}
current['from'] = 'Y'
current['to'] = 'Z'
current['effect'] = YtoZ
addOn.append(current)

current = {}
current['from'] = 'Z'
current['to'] = 'X'
current['effect'] = ZtoX
addOn.append(current)

current = {}
current['from'] = 'Z'
current['to'] = 'Y'
current['effect'] = ZtoY
addOn.append(current)

current = {}
current['from'] = 'Z'
current['to'] = 'Z'
current['effect'] = ZtoZ
addOn.append(current)
inputObject.append(addOn)

#Creates the master input file
fileName = 'masterNodeInputFile'+str(len(inputObject)-2)+'.json'
writeToJson(inputObject,fileName)