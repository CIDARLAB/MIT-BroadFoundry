
def maine() :
    perNcal(nom, count, seq)
    
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

if __name__ == "__maine__" :
    maine()
