from random import random

def maine() :
    reb1_site_chooser(strength, parameterD)
    rap1_site_chooser(strength, parameterD)
    gcr1_site_chooser(strength, parameterD)
    abf1_site_chooser(strength, parameterD)
    mcm1_site_chooser(strength, parameterD)
    rsc3_site(strength, parameterD)
    
def reb1_site_chooser(strength, parameterD):
    
    tf_name = 'REB1'

    reb1_choose = random()
    choice_list = parameterD['reb1_site_choice'][strength]
    
    index = 0
    if reb1_choose <= choice_list[0] :
        index = 1
    if choice_list[0] < reb1_choose <= choice_list[1] :
        index = 2
    if choice_list[1] < reb1_choose <= choice_list[2] :
        index = 3

    reb1_choice = parameterD['reb1'][index]

    return [tf_name, index+1, reb1_choice]

def rap1_site_chooser(strength, parameterD):

    tf_name = 'RAP1'

    rap1_choose = random()
    choice_list = parameterD['rap1_site_choice'][strength]
    
    index = 0
    if rap1_choose <= choice_list[0] :
        index = 1
    if choice_list[0] < rap1_choose <= choice_list[1] :
        index = 2

    rap1_choice = parameterD['rap1'][index]

    return [tf_name, index+1, rap1_choice]

def gcr1_site_chooser(strength, parameterD):
    
    tf_name = 'GCR1'
    
    gcr1_choose = random()
    choice_list = parameterD['gcr1_site_choice'][strength]
    
    index = 0
    if gcr1_choose <= choice_list[0] :
        index = 1
    if choice_list[0] < gcr1_choose <= choice_list[1] :
        index = 2
    if choice_list[1] < gcr1_choose <= choice_list[2] :
        index = 3

    gcr1_choice = parameterD['gcr1'][index]

    return [tf_name, index+1, gcr1_choice]

def abf1_site_chooser(strength, parameterD):
    
    tf_name = 'ABF1'
    
    abf1_choose = random()
    choice_list = parameterD['abf1_site_choice']
    
    index = 0
    if abf1_choose <= choice_list[0] :
        index = 1
    if choice_list[0] < abf1_choose <= choice_list[1] :
        index = 2
        
    abf1_choice = parameterD['abf1'][index]
    
    return [tf_name, index+1, abf1_choice]

def mcm1_site_chooser(strength, parameterD):
    
    tf_name = 'MCM1'

    mcm1_choose = random()
    choice_list = parameterD['mcm1_site_choice']
    
    index = 0
    if mcm1_choose <= choice_list[0] :
        index = 1
        
    mcm1_choice = parameterD['mcm1'][index]
    
    return [tf_name, index+1, mcm1_choice]

def rsc3_site(strength, parameterD):
    
    tf_name  = 'RSC3'
    
    rsc3 = parameterD['rsc3']
    
    return [tf_name, 1, rsc3]

if __name__ == "__maine__" :
    maine()
