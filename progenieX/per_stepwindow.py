from nucleotide_pct_cal import perNcal

def maine() :
    per_stepwindow(input_file, nucpct_file)

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

        step_output.close()


if __name__ == "__maine__" :
    maine()
