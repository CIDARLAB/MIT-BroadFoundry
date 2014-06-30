from random import random
from xlrd import *
from excel_functions import cell

def maine() :
    tss_upstream_chooser(strength, iB, cS)
    tss_el_chooser(strength, iB, cS)
    
def tss_upstream_chooser(strength, iB, cS) :

    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS upstream has three positively correlated elements and 1 negative
    tss_upstream_list = [cell(iB,cS,'R28'), cell(iB,cS,'R29'), cell(iB,cS,'R30'), cell(iB,cS,'R31')]
        
    tssuD = {'VH' : [cell(iB,cS,'I29'), cell(iB,cS,'I30'), cell(iB,cS,'I31')],
             'H' : [cell(iB,cS,'K29'), cell(iB,cS,'K30'), cell(iB,cS,'K31')],
             'M' : [cell(iB,cS,'M29'), cell(iB,cS,'M30'), cell(iB,cS,'M31')],
             'L' : [cell(iB,cS,'O29'), cell(iB,cS,'O30'), cell(iB,cS,'O31')]}

    x = strength
    tss_up_choose = random()
    index = 0
    if tss_up_choose <= tssuD[x][0] :
        index = 1
    if tssuD[x][0] < tss_up_choose <= tssuD[x][1] :
        index = 2
    if tssuD[x][1] < tss_up_choose <= tssuD[x][2] :
        index = 3

    tss_upstream = tss_upstream_list[index]
             
    return tss_upstream

def tss_el_chooser(strength, iB, cS) :
             
    # TSS mutants derived from consensus sequences in Lubliner et al.
    # TSS elements are all positively correlated, but (0) and (1) are more strongly correlated.
    # Therefore I decreased the likelihood of choosing those elements as strength decreased.
    tss_el_list = [cell(iB,cS,'R33'), cell(iB,cS,'R34'), cell(iB,cS,'R35'), cell(iB,cS,'R36')]
             
    tss_elD = {'VH' : [cell(iB,cS,'I34'), cell(iB,cS,'I35'), cell(iB,cS,'I36')],
             'H' : [cell(iB,cS,'K34'), cell(iB,cS,'K35'), cell(iB,cS,'K36')],
             'M' : [cell(iB,cS,'M34'), cell(iB,cS,'M35'), cell(iB,cS,'M36')],
             'L' : [cell(iB,cS,'O34'), cell(iB,cS,'O35'), cell(iB,cS,'O36')]}

    x = strength
    tss_el_choose = random()
    index = 0
    if tss_el_choose <= tss_elD[x][0] :
        index = 1
    if tss_elD[x][0] < tss_el_choose <= tss_elD[x][1] :
        index = 2
    if tss_elD[x][1] < tss_el_choose <= tss_elD[x][2] :
        index = 3

    tss_el = tss_el_list[index]

    return tss_el

if __name__ == "__maine__" :
    maine()
