# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:59:24 2015

@author: Arinze
"""

import json
import inputs

def writeJsonFile(time_axis_params,input_and_logic_gate_dictionaries,fileLoc):
    """
    Writes what was used to a JSON file using pretty print. Saves all
    information needed to reproduce the circuit. Does the opposite process 
    of getFromJsonFile(fileLoc)
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
  
    #data is a JSON array of JSON objects
    #first item is time-axis params (xmin, xmax, #iterations)
    #other elements are dictionaries for genes specifying all properties.
    #supposed to be inputs, repressors, outputs.
    data = json.load(myFile)

    time_axis_params = data[0]
    
    #see also Gate.py:formatJson (could have taken dictionaries directly)
    input_and_logic_gate_dictionaries = data[1:]

    #order is same as dictionary, order matters!    
    input_and_logic_gate_names = []

    
    logic_gate_names = []
    
    for item in input_and_logic_gate_dictionaries:
        input_and_logic_gate_names.append(item['NAME'])
        if item['TYPE']!='INPUT':
            logic_gate_names.append(item['NAME'])
        if item['TYPE']=='INPUT':
            
            #for example: sinInput(t, per=4, amp=1, dis=0, bas=0): ['INPUT'][0] is the function name, [1] is per, [2] is amp, etc.
            #             item['INPUT'] = [inputs.sinInput,4,1,0,0]
            #convert function name to actual function
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