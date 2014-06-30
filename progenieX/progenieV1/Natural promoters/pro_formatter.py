def maine() :
    promoter_formatter()
    
def promoter_formatter():
    anslist = []

    with open ('nat_largeset_promoters.txt') as fasta_input :
        falist = fasta_input.readlines()

        for num, line in enumerate(falist) :
            line = line.replace('\n', '')
            falist[num] = line
            
        seqlist = [line for line in falist if '>' not in line]
        nomlist = [line for line in falist if '>' in line]


        for num, line in enumerate(nomlist) :
            nomlist[num] = line.replace(' ', '_')
            
            ind = num*11
            seq = seqlist[ind]+seqlist[ind+1]+seqlist[ind+2]+seqlist[ind+3]+seqlist[ind+4]+seqlist[ind+5]+seqlist[ind+6]+seqlist[ind+7]+seqlist[ind+8]+seqlist[ind+9]+seqlist[ind+10]
            
            print '%(nom)s\n%(seq)s' % {'nom': nomlist[num], 'seq':seq}
            
            
            

            
              
if __name__ == "__maine__" :
    maine()
