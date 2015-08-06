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
hlX = 10			# half-life (minutes)
aX = np.log(2)/hlX	# degradation rate
B0X = 0.6 			# basal production rate
BX = 600-B0X		# maximum production rate
Xmax = (BX+B0X)/aX 	# maximum x production

hlY = hlX
aY = aX
BY = BX
B0Y = B0X
Ymax = Xmax

n=2 				# hill coefficient
Km=40				# coefficient for activation or repression

inputD = {}
inputD['X']['hl'] = hlX
inputD['X']['a'] = aX
inputD['X']['B0'] = B0X
inputD['X']['B'] = BX
inputD['X']['max'] = Xmax
inputD['X']['n'] = n
inputD['X']['Km'] = Km

inputD['Y']['hl'] = hlY
inputD['Y']['a'] = aY
inputD['Y']['B0'] = B0Y
inputD['Y']['B'] = BY
inputD['Y']['max'] = Ymax
inputD['Y']['n'] = n
inputD['Y']['Km'] = Km
inputObject = [inputD]

# values of -1 = repression, 0 = independent, 1 = activation for autoX/autoY/XtoY/YtoX
# initial conditions set at 0, Xmax/2 and Xmax = B/a
for autoX in range(-1,2):
	for autoY in range(-1,2):
		for XtoY in range(-1,2):
			for YtoX in range(-1,2):
				for iniX in range(3):
					inicX = iniX*Xmax/2.
					for iniY in range(3):
						inicY = iniY*Ymax/2
						current = {}
						current['autoX'] = autoX
						current['autoY'] = autoY
						current['XtoY'] = XtoY
						current['YtoX'] = YtoX
						current['inic'] = [inicX,inicY]
						inputObject.append(current)
						'''
						For individual files
						fileName = 'input_'+str(autoX)+'_'+str(autoY)+'_'+str(XtoY)+'_'+str(YtoX)+'_'+str(iniX)+'_'+str(iniY)+'.json'
						inputObject.writeToJson(fileName)
						'''
#Creates the master input file
fileName = 'masterXYInputFile.json'
writeToJson(inputObject,fileName)
	
