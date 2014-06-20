'''
This set of functions performs common tasks among all scripts.  This includes clearing and opening output files, reading
files of sequences, and outputting properly FASTA formatted sequences.
'''  

def maine() :
    clear(filename)
    lister(FASTA_file)
    fasta_out(seqlist, prefix)
    fastaD_out(seqD, output_file)
    promoter_analysis(prefix)
    
def clear(filename) :

    f = open(filename, 'w')
    f.write('')
    f.close()
    f = open(filename, 'a')

    return f

def lister(FASTA_file) :

    with open(FASTA_file, 'r') as f :
        
        flist = f.readlines()
        
        for n, line in enumerate(flist) :
            
            flist[n] = line.replace('\n', '')

    return flist

def fasta_out(seqlist, prefix) :

    output_filename = '%(prefix)sgen.txt' % {'prefix': prefix}

    output_file = clear(output_filename)

    for num, seq in enumerate(seqlist) :
        output_file.write(">EMY%(prefix)s_%(name)-s\n%(seq)-s\n" % {"prefix": prefix, "name": num+1, "seq": seq})

    output_file.close()
    
def fastaD_out(seqD, prefix) :

    strengths = ['VH', 'H', 'M', 'L']

    output_filename = '%(prefix)sgen.txt' % {'prefix': prefix}

    output_file = clear(output_filename)
    
    for x in strengths :
        
        for num, seq in enumerate(seqD[x]) :
            output_file.write(">EMY%(prefix)s_%(str)-s%(name)-s\n%(seq)-s\n" % {"prefix": prefix,"str": x, "name": num+1, "seq": seqD[x][num]})

    output_file.close()


if __name__ == "__maine__" :
    maine()
