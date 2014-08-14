def maine() :
    clear(filename)
    lister(FASTA_file)
    fasta_out(seqlist, prefix)
    fastaD_out(seqD, output_file)
    
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
    multiplierD = {'VH': 0, 'H': 1, 'M': 2, 'L': 3}
    n_in_strength = len(seqD['VH'])

    output_filename = '%(prefix)sgen.txt' % {'prefix': prefix}

    output_file = clear(output_filename)
    
    for x in strengths :
        
        for num, seq in enumerate(seqD[x]) :
            output_file.write(">EMY%(prefix)s_%(name)-s\n%(seq)-s\n" % {"prefix": prefix,
                                                                        "str": x,
                                                                        "name": num+1+multiplierD[x]*n_in_strength,
                                                                        "seq": seqD[x][num]})

    output_file.close()


if __name__ == "__maine__" :
    maine()
