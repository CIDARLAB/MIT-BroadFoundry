# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:59:24 2015

@author: Arinze
"""

import json
import inputs

def writeJsonFile(time_axis_params,input_and_logic_gate_dictionaries,fileLoc):
    """
    Writes what was used to a JSON file using pretty print. Saves all information
    needed to reporoduce the circuit
    """
    newList = [time_axis_params]
    newList = newList + input_and_logic_gate_dictionaries
    myFile = open(fileLoc,'w')
    json.dump(newList, myFile, sort_keys=True, indent=4, separators=(',', ': '))
    myFile.close()

def getFromJsonFile(fileLoc):
    '''
    Take in a file location and return the values necessary to generate a graph
    '''
    myFile = open(fileLoc,'r')
    data = json.load(myFile)
    time_axis_params = data[0]
    input_and_logic_gate_dictionaries = data[1:]
    input_and_logic_gate_names = []
    logic_gate_names = []
    for item in input_and_logic_gate_dictionaries:
        input_and_logic_gate_names.append(item['NAME'])
        if item['TYPE']!='INPUT':
            logic_gate_names.append(item['NAME'])
        if item['TYPE']=='INPUT':
            x = item['INPUT'][0]
            if x =='inputs.sinInput' or x == 'sinInput':
                item['INPUT'][0] = inputs.sinInput
            elif x =='inputs.linInput' or x == 'linInput':
                item['INPUT'][0] = inputs.linInput
            elif x =='inputs.stepFunction' or x == 'stepFunction':
                item['INPUT'][0] = inputs.stepFunction
            elif x =='inputs.stepInput' or x == 'stepInput':
                item['INPUT'][0] = inputs.stepInput
            elif x =='inputs.squInput' or x == 'squInput':
                item['INPUT'][0] = inputs.squInput
            elif x =='inputs.sawInput' or x == 'sawInput':
                item['INPUT'][0] = inputs.sawInput
            elif x =='inputs.specInput' or x == 'specInput':
                item['INPUT'][0] = inputs.specInput
            else:
                print "This is not a recognized input function."
                return
    return time_axis_params,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names