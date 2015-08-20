# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:00:31 2015

@author: Alex Lim
File that is used to create master input file needed for XY input plots, but can be used 
as basis for converting objects to json files.
"""

import json
import numpy as np

def writeToJson(inputObj,fileLoc):
    """
    Writes the object, preferably composed of lists or dictionaries, to a json file with pretty print
    """
    myFile = open(fileLoc,'w')
    json.dump(inputObj, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()

#Initialization of variables for genes X and Y
hlX = 10		# half-life (minutes)
aX = np.log(2)/hlX	# degradation rate
B0X = 0.6			# basal production rate
BX = 300-B0X		# induced production rate
nX = 2 				# hill coefficient
KmX = 40            # coefficient for activation or repression

hlY = 10
aY = np.log(2)/hlY
B0Y = 0.6
BY = 300-B0Y
nY = 2
KmY = 40

Xmax = (BX+BY+B0X)/aX 	# maximum x production
Ymax = (BX+BY+B0Y)/aY


iniCond = {'X':10.0**-8, 'Y':10.0**3}  
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