'''
This script is part of the promoter analysis group of functions and uses regular expressions to search generated sequences for the motifs that have been substituted or 
arisen by chance.  This data is output into a text file that contains the site information for all sequences.
'''

import re

def maine() :
    tfbs_finder(input_file, tfbs_file)
    reb1_cal(nom,seq)
    rap1_cal(nom,seq)
    gcr1_cal(nom,seq)
    abf1_cal(nom,seq)
    mcm1_cal(nom,seq)
    rsc3_cal(nom,seq)
    pdW_cal(nom,seq)
    
def tfbs_finder(input_file, tfbs_file):
    
    anslist = []

    with open('core_tfbsdata.txt', 'w') as output :
        output.write('')

    with open (input_file) as fasta_input :
        falist = fasta_input.readlines()

        for num, line in enumerate(falist) :
            line = line.replace('\n', '')
            falist[num] = line

        seqlist = [line for line in falist if '>' not in line]
             
        nomlist = [line for line in falist if '>' in line]

        for num, line in enumerate(seqlist) :
            anslist.append(reb1_cal(nomlist[num], line)+rap1_cal(nomlist[num], line)+gcr1_cal(nomlist[num], line)+abf1_cal(nomlist[num],line)+mcm1_cal(nomlist[num], line)+rsc3_cal(nomlist[num], line)+pdW_cal(nomlist[num], line))

        output = open('core_tfbsdata.txt', 'a')
            
        for num, line in enumerate(anslist) :
            ans = anslist[num]
            seqname = ans[0]
            reby1 = ans[1]
            rebl1 = ans[2]
            reby2 = ans[3]
            rebl2 = ans[4]
            reby3 = ans[5]
            rebl3 = ans[6]
            reby4 = ans[7]
            rebl4 = ans[8]
            rapy1 = ans[9]
            rapl1 = ans[10]
            rapy2 = ans[11]
            rapl2 = ans[12]
            rapy3 = ans[13]
            rapl3 = ans[14]
            gcry1 = ans[15]
            gcrl1 = ans[16]
            gcry2 = ans[17]
            gcrl2 = ans[18]
            gcry3 = ans[19]
            gcrl3 = ans[20]
            gcry4 = ans[21]
            gcrl4 = ans[22]
            abfy1 = ans[23]
            abfl1 = ans[24]
            abfy2 = ans[25]
            abfl2 = ans[26]
            abfy3 = ans[27]
            abfl3 = ans[28]
            mcmy1 = ans[29]
            mcml1 = ans[30]
            mcmy2 = ans[31]
            mcml2 = ans[32]
            rscy1 = ans[33]
            rscl1 = ans[34]
            pdWy1 = ans[35]
            pdWl1 = ans[36]
            pdWy2 = ans[37]
            pdWl2 = ans[38]
            pdWy3 = ans[39]
            pdWl3 = ans[40]
            pdWy4 = ans[41]
            pdWl4 = ans[42]
            output.write("%(N)s %(reby1)s %(rebl1)s %(reby2)s %(rebl2)s %(reby3)s %(rebl3)s %(reby4)s %(rebl4)s \
%(rapy1)s %(rapl1)s %(rapy2)s %(rapl2)s %(rapy3)s %(rapl3)s \
%(gcry1)s %(gcrl1)s %(gcry2)s %(gcrl2)s %(gcry3)s %(gcrl3)s \
%(gcry4)s %(gcrl4)s %(abfy1)s %(abfl1)s %(abfy2)s %(abfl2)s \
%(abfy3)s %(abfl3)s %(mcmy1)s %(mcml1)s %(mcmy2)s %(mcml2)s \
%(rscy1)s %(rscl1)s %(pdWy1)s %(pdWl1)s %(pdWy2)s %(pdWl2)s \
%(pdWy3)s %(pdWl3)s %(pdWy4)s %(pdWl4)s\n" \
            % {'N':seqname, 'reby1':reby1, 'rebl1': rebl1, 'reby2':reby2, 'rebl2': rebl2, 'reby3':reby3, \
               'rebl3': rebl3, 'reby4':reby4, 'rebl4': rebl4, 'rapy1':rapy1, 'rapl1': rapl1, 'rapy2':rapy2, \
               'rapl2': rapl2, 'rapy3':rapy3, 'rapl3': rapl3, 'gcry1':gcry1, 'gcrl1': gcrl1, 'gcry2':gcry2, \
               'gcrl2': gcrl2, 'gcry3':gcry3, 'gcrl3': gcrl3, 'gcry4': gcry4, 'gcrl4': gcrl4, 'abfy1':abfy1, \
               'abfl1': abfl1, 'abfy2':abfy2, 'abfl2': abfl2, 'abfy3': abfy3, 'abfl3': abfl3, 'mcmy1':mcmy1,\
               'mcml1': mcml1, 'mcmy2':mcmy2, 'mcml2': mcml2, 'rscy1': rscy1, 'rscl1': rscl1, 'pdWy1':pdWy1, \
               'pdWl1': pdWl1, 'pdWy2':pdWy2, 'pdWl2': pdWl2, 'pdWy3':pdWy3, 'pdWl3': pdWl3, 'pdWy4':pdWy4, \
               'pdWl4': pdWl4})
        output.close()

def reb1_cal(nom, seq) :

    y = [0,0,0,0]
    l = [0,0,0,0]
    m = [0,0,0,0]
    sites = [0,1,2,3]
    siteD = {0 : '(TTACCCGT)', 1 : '(TCACCCGT)', 2 : '(CAGCCCTT)', 3 : '(TTACCCGG)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return nom, y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]

def rap1_cal(nom, seq) :

    y = [0,0,0]
    l = [0,0,0]
    m = [0,0,0]
    sites = [0,1,2]
    siteD = {0 : '(ACACCCAAGCAT)', 1 : '(ACACCTGGACAT)', 2 : '(ACCCCTTTTTTAC)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2]

def gcr1_cal(nom, seq) :

    y = [0,0,0,0]
    l = [0,0,0,0]
    m = [0,0,0,0]
    sites = [0,1,2,3]
    siteD = {0 : '(CGACTTCCT)', 1 : '(CGGCATCCA)', 2 : '(CAGCTTCCT)', 3 : '(CAACGGAAG)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]

def abf1_cal(nom, seq) :

    y = [0,0,0]
    l = [0,0,0]
    m = [0,0,0]
    sites = [0,1,2]
    siteD = {0 : '(AGCCGTAAATAGTTATCTTCCAAG)', 1 : '(ATCATCTATCACG)', 2 : '(GTCATTTTACACG)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return nom, y[0], l[0], y[1], l[1], y[2], l[2]

def mcm1_cal(nom, seq) :

    y = [0,0]
    l = [0,0]
    m = [0,0]
    sites = [0,1]
    siteD = {0 : '(TTTCCGAAAACGGAAAT)', 1 : '(ATACCAAATACGGTAAT)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1]

def rsc3_cal(nom, seq) :

    y = [0]
    l = [0]
    m = [0]
    sites = [0]
    siteD = {0 : '(CGCGC)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0]

def pdW_cal(nom, seq) :

    y = [0,0,0,0]
    l = [0,0,0,0]
    m = [0,0,0,0]
    sites = [0,1,2,3]
    siteD = {0 : '(TTTTTTTTTTTTT)', 1 : '(AAAAAAAAAAAAA)', 2 : '(TTAATTTAATTTT)', 3:'([AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT])'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]
   
if __name__ == "__maine__" :
    maine()
