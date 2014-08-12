from random import random
from xlrd import *
from excel_functions import cell

def maine() :
    reb1_site_chooser(strength, parameterD)
    rap1_site_chooser(strength, parameterD)
    gcr1_site_chooser(strength, parameterD)
    abf1_site_chooser(strength, parameterD)
    mcm1_site_chooser(strength, parameterD)
    rsc3_site(strength, parameterD)
    
def reb1_site_chooser(strength, parameterD):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    tf_name = 'REB1'
    
    # BINDING SITES FROM YEASTRACT AND MOGNO ET AL.
    reb1 = [cell(iB,uS,'S5'), cell(iB,uS,'S6'), cell(iB,uS,'S7'), cell(iB,uS,'S8')]

    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # for VH strength, only the first or second REB1 site is possible to choose.  For H strength, three are possible, adding
    # in the consensus site published in YEASTRACT.  For M, the 'weak' site from Mogno is included as a choice. For L, the
    # weak site is most likely.
    reb1_site_choiceD = {'VH' : [cell(iB,uS,'J6'), cell(iB,uS,'J7'), cell(iB,uS,'J8')],
                         'H' : [cell(iB,uS,'L6'), cell(iB,uS,'L7'), cell(iB,uS,'L8')],
                         'M' : [cell(iB,uS,'N6'), cell(iB,uS,'N7'), cell(iB,uS,'N8')],
                         'L' : [cell(iB,uS,'P6'), cell(iB,uS,'P7'), cell(iB,uS,'P8')]}

    reb1_choose = random()
    index = 0
    if reb1_choose <= reb1_site_choiceD[strength][0] :
        index = 1
    if reb1_site_choiceD[strength][0] < reb1_choose <= reb1_site_choiceD[strength][1] :
        index = 2
    if reb1_site_choiceD[strength][1] < reb1_choose <= reb1_site_choiceD[strength][2] :
        index = 3

    reb1_choice = reb1[index]

    return [tf_name, index+1, reb1_choice]

def rap1_site_chooser(strength, parameterD):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    tf_name = 'RAP1'
    
    # BINDING SITES FROM YEASTRACT AND MOGNO ET AL.
    rap1 = [cell(iB,uS,'S14'), cell(iB,uS,'S15'), cell(iB,uS,'S16')]

    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # Very similar to the code for the REB1_chooser, so check that out for more discussion.
    rap1_site_choiceD = {'VH' : [cell(iB,uS,'J15'), cell(iB,uS,'J16')],
                         'H' : [cell(iB,uS,'L15'), cell(iB,uS,'L16')],
                         'M' : [cell(iB,uS,'N15'), cell(iB,uS,'N16')],
                         'L' : [cell(iB,uS,'P15'), cell(iB,uS,'P16')]}

    rap1_choose = random()
    index = 0
    if rap1_choose <= rap1_site_choiceD[strength][0] :
        index = 1
    if rap1_site_choiceD[strength][0] < rap1_choose <= rap1_site_choiceD[strength][1] :
        index = 2

    rap1_choice = rap1[index]

    return [tf_name, index+1, rap1_choice]

def gcr1_site_chooser(strength, parameterD):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    tf_name = 'GCR1'
    
    # BINDING SITES FROM YEASTRACT, MOGNO ET AL., AND BAI ET AL.
    gcr1 = [cell(iB,uS,'S22'), cell(iB,uS,'S23'), cell(iB,uS,'S24'), cell(iB,uS,'S25')]
    
    # This dictionary changes the probability of picking a strong or weak binding site depending on the strength of the UAS.
    # Very similar to the code for the REB1_chooser, so check that out for more discussion.
    gcr1_site_choiceD = {'VH' : [cell(iB,uS,'J23'), cell(iB,uS,'J24'), cell(iB,uS,'J25')],
                         'H' : [cell(iB,uS,'L23'), cell(iB,uS,'L24'), cell(iB,uS,'L25')],
                         'M' : [cell(iB,uS,'N23'), cell(iB,uS,'N24'), cell(iB,uS,'N25')],
                         'L' : [cell(iB,uS,'P23'), cell(iB,uS,'P24'), cell(iB,uS,'P25')]}

    gcr1_choose = random()
    index = 0
    if gcr1_choose <= gcr1_site_choiceD[strength][0] :
        index = 1
    if gcr1_site_choiceD[strength][0] < gcr1_choose <= gcr1_site_choiceD[strength][1] :
        index = 2
    if gcr1_site_choiceD[strength][1] < gcr1_choose <= gcr1_site_choiceD[strength][2] :
        index = 3

    gcr1_choice = gcr1[index]

    return [tf_name, index+1, gcr1_choice]

def abf1_site_chooser(strength, parameterD):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    tf_name = 'ABF1'
    
    # BINDING SITES FROM YEASTRACT, LAST TWO BASED OFF "RTCRYYYNNNACG" 
    abf1 = [cell(iB,uS,'S30'), cell(iB,uS,'S31'), cell(iB,uS,'S32')]

    abf1_choose = random()
    index = 0
    if abf1_choose <= cell(iB,uS,'J31') :
        index = 1
    if 0.33 < abf1_choose <= cell(iB,uS,'J32') :
        index = 2
        
    abf1_choice = abf1[index]
    
    return [tf_name, index+1, abf1_choice]

def mcm1_site_chooser(strength, parameterD):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    tf_name = 'MCM1'
    
    # BINDING SITES FROM BAI ET AL.
    mcm1 = [cell(iB,uS,'S37'), cell(iB,uS,'S38')]
    
    mcm1_choose = random()
    index = 0
    if mcm1_choose <= cell(iB,uS,'J38') :
        index = 1
        
    mcm1_choice = mcm1[index]
    
    return [tf_name, index+1, mcm1_choice]

def rsc3_site(strength, parameterD):
    # Define location of data 
    iB = 'ProGenie_Parameters.xlsx'
    gS = 'General'
    cS = 'Core'
    uS = 'UAS'
    
    tf_name  = 'RSC3'
    
    # BINDING SITE FROM BAI ET AL.
    rsc3 = cell(iB,uS,'S41')
    
    return [tf_name, 1, rsc3]

if __name__ == "__maine__" :
    maine()
