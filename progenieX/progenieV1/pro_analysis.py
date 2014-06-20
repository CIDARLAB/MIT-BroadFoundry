'''
This function is a simple way of activating all functions that analyze sequences created by PRO_genie.
'''

from per_stepwindow import per_stepwindow
from tata_stepwindow import tata_stepwindow
from tfbs_finder import tfbs_finder
from tss_finder import tss_finder

def maine() :
    promoter_analysis(prefix)
    
def promoter_analysis(prefix) :

    input_file = '%(prefix)sgen.txt' % {'prefix': prefix}
    nucpct_file = '%(prefix)s_stepdata.txt' % {'prefix': prefix}
    tata_file = '%(prefix)s_TATAdata.txt' % {'prefix': prefix}
    tfbs_file = '%(prefix)s_tfbsdata.txt' % {'prefix': prefix}
    tss_file = '%(prefix)s_tssdata.txt' % {'prefix': prefix}

    per_stepwindow(input_file, nucpct_file)
    tata_stepwindow(input_file, tata_file)
    tfbs_finder(input_file, tfbs_file)
    tss_finder(input_file, tss_file)


if __name__ == "__maine__" :
    maine()
