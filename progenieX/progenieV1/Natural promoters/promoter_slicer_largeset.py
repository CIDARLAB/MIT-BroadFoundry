def maine() :
    promoter_slicer()
    
def promoter_slicer():

    typL = ['uas2', 'uas1', 'core']
    sliceD = {'uas2': [450,300], 'uas1': [300,150], 'core': [150, 0]}
    strengths = ['VH', 'H', 'M', 'L']
    strD = {'VH': [0,22], 'H': [23,49], 'M': [50,78], 'L': [79,124]}

    anslist = []
        
    A = 'GTGC'
    B = 'AATG'
    J1 = 'TTCT'
    J2 = 'AAAC'
    J3 = 'ACTA'
    J4 = 'CTTA'
    J5 = 'CCGA'
    J6 = 'GATA'
    J7 = 'CCGA'
    J8 = 'GACC'

    scarD = {'VH' : {'uas2':[A, J1], 'uas1':[J1,J2], 'core':[J2,B]},
             'H' : {'uas2':[A, J3], 'uas1':[J3,J4], 'core':[J4,B]},
             'M' : {'uas2':[A, J5], 'uas1':[J5,J6], 'core':[J6,B]},
             'L' : {'uas2':[A, J7], 'uas1':[J7,J8], 'core':[J8,B]}}

    for y in typL :
        output = open('%(typ)s_lg_natpros.txt' % {'typ':y}, 'w')
        output.write('')
        output.close()
    
    with open ('nat_largeset_promoters.txt') as fasta_input :

        falist = fasta_input.readlines()

        for num, line in enumerate(falist) :
            line = line.replace('\n', '')
            falist[num] = line
            
        seqlist = [line for line in falist if '>' not in line]
        nomlist = [line for line in falist if '>' in line]

        for x in strengths :
            
            start_i = strD[x][0]
            end_i = strD[x][1]
            
            while start_i <= end_i :

                for y in typL :

                    promoter = seqlist[start_i]
                        
                    Lscar, Rscar = scarD[x][y][0], scarD[x][y][1]
            
                    seq = Lscar+promoter[len(promoter)-sliceD[y][0]:len(promoter)-sliceD[y][1]]+Rscar

                    output = open('%(typ)s_lg_natpros.txt' % {'typ':y}, 'a')
            
                    output.write('%(n)s_%(typ)s_%(str)s\n%(seq)s\n' % {'n':nomlist[start_i],'typ':y,'str':x,'seq':seq})

                    output.close()

                start_i = start_i+1
                
            
              
if __name__ == "__maine__" :
    maine()
