import re

def maine() :
    reb1_cal(nom,seq)
    rap1_cal(nom,seq)
    gcr1_cal(nom,seq)
    abf1_cal(nom,seq)
    mcm1_cal(nom,seq)
    rsc3_cal(nom,seq)
    pdW_cal(nom,seq)

def reb1_cal(nom, seq) :

    y = [0,0,0,0]
    l = [0,0,0,0]
    m = [0,0,0,0]
    sites = [0,1,2,3]
    siteD = {0 : '(TTACCCGT)', 1 : '(TCACCCGT)', 2 : '(CAGCCCTT)', 3 : '(TTACCCGG)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return nom, y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]

def rap1_cal(nom, seq) :

    y = [0,0,0]
    l = [0,0,0]
    m = [0,0,0]
    sites = [0,1,2]
    siteD = {0 : '(ACACCCAAGCAT)', 1 : '(ACACCTGGACAT)', 2 : '(ACCCCTTTTTTAC)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2]

def gcr1_cal(nom, seq) :

    y = [0,0,0,0]
    l = [0,0,0,0]
    m = [0,0,0,0]
    sites = [0,1,2,3]
    siteD = {0 : '(CGACTTCCT)', 1 : '(CGGCATCCA)', 2 : '(CAGCTTCCT)', 3 : '(CAACGGAAG)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]

def abf1_cal(nom, seq) :

    y = [0,0,0]
    l = [0,0,0]
    m = [0,0,0]
    sites = [0,1,2]
    siteD = {0 : '(AGCCGTAAATAGTTATCTTCCAAG)', 1 : '(ATCATCTATCACG)', 2 : '(GTCATTTTACACG)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return nom, y[0], l[0], y[1], l[1], y[2], l[2]

def mcm1_cal(nom, seq) :

    y = [0,0]
    l = [0,0]
    m = [0,0]
    sites = [0,1]
    siteD = {0 : '(TTTCCGAAAACGGAAAT)', 1 : '(ATACCAAATACGGTAAT)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1]

def rsc3_cal(nom, seq) :

    y = [0]
    l = [0]
    m = [0]
    sites = [0]
    siteD = {0 : '(CGCGC)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0]

def pdW_cal(nom, seq) :

    y = [0,0,0,0]
    l = [0,0,0,0]
    m = [0,0,0,0]
    sites = [0,1,2,3]
    siteD = {0 : '(TTTTTTTTTTTTT)', 1 : '(AAAAAAAAAAAAA)', 2 : '(TTAATTTAATTTT)', 3:'([AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT][AT])'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3]
   
if __name__ == "__maine__" :
    maine()
