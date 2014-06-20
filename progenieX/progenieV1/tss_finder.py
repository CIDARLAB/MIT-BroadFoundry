'''
This function searches for all possible combinations of the transcription start site elements that are possible in
PRO_genie.  It is part of the promoter analysis group of functions.
'''

import re

def maine() :
    tss_finder(input_file, tss_file)
    tss_cal(nom,seq)

def tss_finder(input_file, tss_file):
    
    anslist = []

    with open(tss_file, 'w') as output :
        output.write('')

    with open (input_file) as fasta_input :
        falist = fasta_input.readlines()

        for num, line in enumerate(falist) :
            line = line.replace('\n', '')
            falist[num] = line

        seqlist = [line for line in falist if '>' not in line]
             
        nomlist = [line for line in falist if '>' in line]

        for num, line in enumerate(seqlist) :
            anslist.append(tss_cal(nomlist[num], line))

        output = open(tss_file, 'a')
            
        for num, line in enumerate(anslist) :
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
               'tssl15': tssl15, 'tssy16':tssy16, 'tssl16': tssl16})

        output.close()
        
def tss_cal(nom, seq) :

    y = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    l = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    m = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    sites = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    siteD = {0 : '(TTTTCAAA)', 1 : '(TTTTCAAT)', 2 : '(TTTTCACC)', 3 : '(TTTTACAA)', \
             4 : '(TTCTCAAA)', 5 : '(TTCTCAAT)', 6 : '(TTCTCACC)', 7 : '(TTCTACAA)', \
             8 : '(CTTACAAA)', 9 : '(CTTACAAT)', 10 : '(CTTACACC)', 11 : '(CTTAACAA)', \
             12 : '(AGCGCAAA)', 13 : '(AGCGCAAT)', 14 : '(AGCGCACC)', 15 : '(AGCGACAA)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return nom, y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3], y[4], l[4], y[5], l[5], y[6], l[6], y[7], l[7], y[8], l[8], y[9], l[9], y[10], l[10], y[11], l[11], y[12], l[12], y[13], l[13], y[14], l[14], y[15], l[15]

  
if __name__ == "__maine__" :
    maine()
