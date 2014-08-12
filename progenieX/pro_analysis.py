import re

def maine():
    promoter_analysis(raw_input('Prefix:'))
    
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
    
def per_stepwindow(input_file, nucpct_file):
    
    anslist = []
    
    with open(nucpct_file, 'w') as output :
        output.write('')
        
    with open (input_file) as fasta_input :
        falist = fasta_input.readlines()

        for num, line in enumerate(falist) :
            line = line.replace('\n', '')
            falist[num] = line

        seqlist = [line for line in falist if '>' not in line]
             
        nomlist = [line for line in falist if '>' in line]

        for num, line in enumerate(seqlist) :
            step = 10
            sliceL = 0
            sliceR = 20
            count = 1
            while sliceR <= len(line) :
                anslist.append(perNcal(nomlist[num], count, line[sliceL:sliceR]))
                #TATA
                '''anslist.append(tatafind(nomlist[num], line[sliceL:sliceR], count, step))'''

                #TFBS
                '''anslist.append(reb1_cal(nomlist[num], line)+
                               rap1_cal(nomlist[num], line)+
                               gcr1_cal(nomlist[num], line)+
                               abf1_cal(nomlist[num],line)+
                               mcm1_cal(nomlist[num], line)+
                               rsc3_cal(nomlist[num], line)+
                               pdW_cal(nomlist[num], line))'''

                #TSS
                '''anslist.append(tss_cal(nomlist[num], line))'''

                sliceL = sliceL + step
                sliceR = sliceR + step
                count = count + 1
                    
        step_output = open(nucpct_file, 'a')
            
        for num, line in enumerate(anslist) :
            ans = anslist[num]
            seqname = ans[0]
            count = ans[1]
            pct = ans[2]
            pctA = pct[0]
            pctT = pct[1]
            pctC = pct[2]
            pctG = pct[3]
            step_output.write("%(N)s %(ct)s %(A)s %(T)s %(C)s %(G)s\n" % {'N':seqname,'ct': count, 'A':pctA,'T':pctT,'C':pctC,'G':pctG})

        #TATA
        '''for num, line in enumerate(anslist) :
            ans = anslist[num]
            seqname = ans[0]
            tataF = ans[1]
            locationF = ans[2]
            tataR = ans[3]
            locationR = ans[4]
            
            output.write("%(N)s,%(tataF)s,%(locF)s,%(tataR)s,%(locR)s\n" % {'N':seqname, 'tataF':tataF, 'locF': locationF, 'tataR':tataR, 'locR': locationR})'''

        #TFBS
        '''for num, line in enumerate(anslist) :
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
               'pdWl4': pdWl4})'''

        #TSS
        '''for num, line in enumerate(anslist) :
            ans = anslist[num]
            seqname = ans[0]
            tssy1 = ans[1]
            tssl1 = ans[2]
            tssy2 = ans[3]
            tssl2 = ans[4]
            tssy3 = ans[5]
            tssl3 = ans[6]
            tssy4 = ans[7]
            tssl4 = ans[8]
            tssy5 = ans[9]
            tssl5 = ans[10]
            tssy6 = ans[11]
            tssl6 = ans[12]
            tssy7 = ans[13]
            tssl7 = ans[14]
            tssy8 = ans[15]
            tssl8 = ans[16]
            tssy9 = ans[17]
            tssl9 = ans[18]
            tssy10 = ans[19]
            tssl10 = ans[20]
            tssy11 = ans[21]
            tssl11 = ans[22]
            tssy12 = ans[23]
            tssl12 = ans[24]
            tssy13 = ans[25]
            tssl13 = ans[26]
            tssy14 = ans[27]
            tssl14 = ans[28]
            tssy15 = ans[29]
            tssl15 = ans[30]
            tssy16 = ans[31]
            tssl16 = ans[32]
            output.write("%(N)s %(tssy1)s %(tssl1)s %(tssy2)s %(tssl2)s %(tssy3)s %(tssl3)s %(tssy4)s %(tssl4)s \
%(tssy5)s %(tssl5)s %(tssy6)s %(tssl6)s %(tssy7)s %(tssl7)s \
%(tssy8)s %(tssl8)s %(tssy9)s %(tssl9)s %(tssy10)s %(tssl10)s \
%(tssy11)s %(tssl11)s %(tssy12)s %(tssl12)s %(tssy13)s %(tssl13)s \
%(tssy14)s %(tssl14)s %(tssy15)s %(tssl15)s %(tssy16)s %(tssl16)s\n" \
            % {'N':seqname, 'tssy1':tssy1, 'tssl1': tssl1, 'tssy2':tssy2, 'tssl2': tssl2, 'tssy3':tssy3, \
               'tssl3': tssl3, 'tssy4':tssy4, 'tssl4': tssl4, 'tssy5':tssy5, 'tssl5': tssl5, 'tssy6':tssy6, \
               'tssl6': tssl6, 'tssy7':tssy7, 'tssl7': tssl7, 'tssy8':tssy8, 'tssl8': tssl8, 'tssy9':tssy9, \
               'tssl9': tssl9, 'tssy10':tssy10, 'tssl10': tssl10, 'tssy11':tssy11, 'tssl11': tssl11, 'tssy12':tssy12, \
               'tssl12': tssl12, 'tssy13':tssy13, 'tssl13': tssl13, 'tssy14':tssy14, 'tssl14': tssl14, 'tssy15':tssy15, \
               'tssl15': tssl15, 'tssy16':tssy16, 'tssl16': tssl16})'''

        step_output.close()

def perNcal(nom, count, seq) :
    
    # Define iterating list ATCG
    ATCG = ["A", "T", "C", "G"]

    # Define dictionaries
    cD = {}
    pD = {}
    outlist = []
    
    # Set all nucleotides to uppercase
    sEQ = seq.upper()

    # Calculate length of the input sequence
    seqlen = float(len(sEQ))
  
    # Iteration loop for calculating the percentage of each base
    for x in ATCG :
        
        cD[x] = sEQ.count(x)

        pD[x] = cD[x] / seqlen * 100

        outlist.append("%5.3f" % (pD[x]))

    return nom, count, outlist

def tatafind(nom, seq, count, step) :

    yesF, yesR = 0, 0
    locationF, locationR = 'NA', 'NA'
    
    match_F = re.search('(TATA[AT]A[AT][AG])', seq)
    match_R = re.search('([CT][AT]T[AT]TATA)', seq)
    
    if match_F:
        yesF = 1
        loc_loc = match_F.start()
        locationF = loc_loc+count*step

    if match_R:
        yesR = 1
        loc_loc = match_R.start()
        locationR = loc_loc+count*step

    return nom, yesF, locationF, yesR, locationR

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
    siteD = {0 : '(TTTTTTTTTTTTT)',
             1 : '(AAAAAAAAAAAAA)',
             2 : '(TTAATTTAATTTT)',
             3:'([AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT])'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]

def tss_cal(nom, seq) :

    y = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    l = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    m = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    
    sites = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    
    siteD = {0 : '(TTTTCAAA)', 1 : '(TTTTCAAT)', 2 : '(TTTTCACC)', 3 : '(TTTTACAA)', 
             4 : '(TTCTCAAA)', 5 : '(TTCTCAAT)', 6 : '(TTCTCACC)', 7 : '(TTCTACAA)', 
             8 : '(CTTACAAA)', 9 : '(CTTACAAT)', 10 : '(CTTACACC)', 11 : '(CTTAACAA)', 
             12 : '(AGCGCAAA)', 13 : '(AGCGCAAT)', 14 : '(AGCGCACC)', 15 : '(AGCGACAA)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return [nom,
            y[0], l[0], y[1], l[1],
            y[2], l[2], y[3], l[3],
            y[4], l[4], y[5], l[5],
            y[6], l[6], y[7], l[7],
            y[8], l[8], y[9], l[9],
            y[10], l[10], y[11], l[11],
            y[12], l[12], y[13], l[13],
            y[14], l[14], y[15], l[15]]

if __name__ == "__maine__" :
    maine()
