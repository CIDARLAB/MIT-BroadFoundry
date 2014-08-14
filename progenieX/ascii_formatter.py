from common_functions import *

def maine() :
    ascii_formatter()
    
def ascii_formatter() :

    coreL = lister('coregen.txt')
    uas1L = lister('uas1gen.txt')
    uas2L = lister('uas2gen.txt')

    cores = [line for line in coreL if '>' not in line]
    uas1s = [line for line in uas1L if '>' not in line]
    uas2s = [line for line in uas2L if '>' not in line]

    print len(cores)
    
    clear('EMY_ascii_08142014.txt')
    
    with open('EMY_ascii_08142014.txt', 'a') as f :

        for n, seq in enumerate(cores) :
            f.write('%(seq)s\n' % {'seq' : seq})
            f.write('%(seq)s\n' % {'seq' : uas1s[n]})
            f.write('%(seq)s\n' % {'seq' : uas2s[n]})


if __name__ == "__maine__" :
    maine()
