# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
File that creates master input file containing initial conditions, parameters 
and vectors for the nodes X and Y.

Uses the method writeToJson to convert a list of dictionary and list objects 
into json format.

Parameters came from the paper: 
A synthetic oscillatory network of transcriptional regulators (Elowitz and 
Leibler, 2000).
"""

import json
import numpy as np


def writeToJson(inputObj,fileLoc):
    """
    Writes the object, preferably a list composed of lists or dictionaries, to a json file with pretty print
    """
    myFile = open(fileLoc,'w')
    json.dump(inputObj, myFile, sort_keys=True, indent=4, separators=(',', ': '))

#Initialization of variables for genes X and Y
hlX = 10		# half-life (minutes)
aX = np.log(2)/hlX	# degradation rate
KmX = 40            # coefficient for activation or repression
B0X = 0.6			# basal production rate
BX = 600-B0X		# maximum production rate
Xmax = (BX+B0X)/aX 	# maximum x production
nX = 2 				# hill coefficient

hlY = hlX
aY = aX
BY = BX
B0Y = B0X
Ymax = Xmax
nY = nX
KmY = 40

iniCond = {'X':0.0, 'Y':0.0}  
inputObject = [iniCond]

inputD = {'X':{},'Y':{}}
inputD['X']['name'] = 'X'
inputD['X']['hl'] = hlX
inputD['X']['a'] = aX
inputD['X']['B0'] = B0X
inputD['X']['B'] = BX
inputD['X']['max'] = Xmax
inputD['X']['n'] = nX
inputD['X']['Km'] = KmX

inputD['Y']['name'] = 'Y'
inputD['Y']['hl'] = hlY
inputD['Y']['a'] = aY
inputD['Y']['B0'] = B0Y
inputD['Y']['B'] = BY
inputD['Y']['max'] = Ymax
inputD['Y']['n'] = nY
inputD['Y']['Km'] = KmY
inputObject.append(inputD)
'''
# Creates vector object of X-Y node interactions (activation, repression and none); size=64 
for XtoX in range(-1,2):
	for XtoY in range(-1,2):
		for YtoX in range(-1,2):
			for YtoY in range(-1,2):
				addOn = []
				current = {}
				current['to'] = 'X'
				current['from'] = 'X'
				current['effect'] = XtoX
				addOn.append(current)

				current = {}
				current['to'] = 'X'
				current['from'] = 'Y'
				current['effect'] = XtoY
				addOn.append(current)

				current = {}
				current['to'] = 'Y'
				current['from'] = 'X'
				current['effect'] = YtoX
				addOn.append(current)

				current = {}
				current['to'] = 'Y'
				current['from'] = 'Y'
				current['effect'] = YtoY
				addOn.append(current)

				inputObject.append(addOn)
'''
'''
# Creates vector object of X-Y node interactions (repression and none); size=16
for XtoX in range(-1,1):
	for XtoY in range(-1,1):
		for YtoX in range(-1,1):
			for YtoY in range(-1,1):
				addOn = []
				current = {}
				current['to'] = 'X'
				current['from'] = 'X'
				current['effect'] = XtoX
				addOn.append(current)

				current = {}
				current['to'] = 'X'
				current['from'] = 'Y'
				current['effect'] = XtoY
				addOn.append(current)

				current = {}
				current['to'] = 'Y'
				current['from'] = 'X'
				current['effect'] = YtoX
				addOn.append(current)

				current = {}
				current['to'] = 'Y'
				current['from'] = 'Y'
				current['effect'] = YtoY
				addOn.append(current)

				inputObject.append(addOn)
'''

# Creates an individual case of X-Y node interactions.
XtoX = 0
XtoY = 1
YtoX = -1
YtoY = 0
addOn = []
current = {}
current['to'] = 'X'
current['from'] = 'X'
current['effect'] = XtoX
addOn.append(current)

current = {}
current['to'] = 'Y'
current['from'] = 'X'
current['effect'] = XtoY
addOn.append(current)

current = {}
current['to'] = 'X'
current['from'] = 'Y'
current['effect'] = YtoX
addOn.append(current)

current = {}
current['to'] = 'Y'
current['from'] = 'Y'
current['effect'] = YtoY
addOn.append(current)

inputObject.append(addOn)


#Creates the master input file
fileName = 'masterXYInputFile'+str(len(inputObject)-2)+'.json'
writeToJson(inputObject,fileName)