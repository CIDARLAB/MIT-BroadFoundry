'''
This function searches for consensus TATA sequences in a stepwise manner through a sequence.  It is part of the 
promoter analysis group of functions.
'''

from tata_find import tatafind

def maine() :
    tata_stepwindow(input_file, tata_file)
    
def tata_stepwindow(input_file, tata_file):
    
    anslist = []

    with open(tata_file, 'w') as output :
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
                    anslist.append(tatafind(nomlist[num], line[sliceL:sliceR], count, step))
                    sliceL = sliceL + step
                    sliceR = sliceR + step
                    count = count + 1

        output = open(tata_file, 'a')
            
        for num, line in enumerate(anslist) :
            ans = anslist[num]
            seqname = ans[0]
            tataF = ans[1]
            locationF = ans[2]
            tataR = ans[3]
            locationR = ans[4]
            
            output.write("%(N)s,%(tataF)s,%(locF)s,%(tataR)s,%(locR)s\n" % {'N':seqname, 'tataF':tataF, 'locF': locationF, 'tataR':tataR, 'locR': locationR})

        output.close()
        
if __name__ == "__maine__" :
    maine()
