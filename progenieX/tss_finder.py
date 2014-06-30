from tss_cal import tss_cal

def maine() :
    tss_finder(input_file, tss_file)

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
   
if __name__ == "__maine__" :
    maine()
