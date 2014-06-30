import re

def maine() :
    tss_cal(nom, seq)
    
def tss_cal(nom, seq) :

    y = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    l = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    m = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    sites = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    siteD = {0 : '(TTTTCAAA)', 1 : '(TTTTCAAT)', 2 : '(TTTTCACC)', 3 : '(TTTTACAA)', \
             4 : '(TTCTCAAA)', 5 : '(TTCTCAAT)', 6 : '(TTCTCACC)', 7 : '(TTCTACAA)', \
             8 : '(CTTACAAA)', 9 : '(CTTACAAT)', 10 : '(CTTACACC)', 11 : '(CTTAACAA)', \
             12 : '(AGCGCAAA)', 13 : '(AGCGCAAT)', 14 : '(AGCGCACC)', 15 : '(AGCGACAA)'}

    for x in sites : 
        y[x] = 0
        l[x] = 'NA'
        m[x] = re.search(siteD[x], seq)
    
        if m[x]:
            y[x] = 1
            l[x] = m[x].start()
        
    return nom, y[0], l[0], y[1], l[1], y[2], l[2], y[3], l[3], y[4], l[4], y[5], l[5], y[6], l[6], y[7], l[7], y[8], l[8], y[9], l[9], y[10], l[10], y[11], l[11], y[12], l[12], y[13], l[13], y[14], l[14], y[15], l[15]

   
if __name__ == "__maine__" :
    maine()
