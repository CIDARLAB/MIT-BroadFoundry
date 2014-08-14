from random import random
from reverser import rev_comp

def maine() :
    re_eraser(raw_input())
    atg_eraser(seq, strength)
    nab_nrd_eraser(seq, strength)

def re_eraser(seq, parameterD) :

    reD = parameterD['re']
    
    re_count = 0
    
    for k in reD:
        
        if reD[k]['seq'] in seq :
            
            seq = seq.replace(reD[k]['seq'],
                              reD[k]['fix'])
            re_count = re_count + 1
            
    return seq

def atg_eraser(seq, strength, parameterD) :

    if strength is not "L":
    
        choose_fix = random()
        choice_list = parameterD['atg_erase']['fix_choice']
        
        index = 0
        if 0 <= choose_fix <= choice_list[0]:
            index = 1
        if choice_list[0] <= choose_fix <= choice_list[0]:
            index = 2
            
        atg = parameterD['atg_erase']['atg']
        fix = parameterD['atg_erase']['fix'][index]
        
        if atg in seq :
            seq = seq.replace(atg, fix)
            
    return [seq, 1]

def nab_nrd_eraser(seq, strength, parameterD) :
    
    if strength is not "L":
        
        if parameterD['nn_erase']['nab3'] in seq :
            seq = seq.replace(parameterD['nn_erase']['nab3'],
                              parameterD['nn_erase']['fix'][0])
            
        if parameterD['nn_erase']['nrd1_1'] in seq :
            seq = seq.replace(parameterD['nn_erase']['nrd1_1'],
                              parameterD['nn_erase']['fix'][1])
            
        if parameterD['nn_erase']['nrd1_2'] in seq :
            seq = seq.replace(parameterD['nn_erase']['nrd1_2'],
                              parameterD['nn_erase']['fix'][1])
            
    return [seq, 1]

if __name__ == "__maine__" :
    maine()
