import Reads
import time
import copy
import matplotlib.pyplot as plt

#File Locations
genomelocation = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/MG1655NC_000913.fasta'
fp_rdm_f = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_fp_rdm_pooled_f.wig'
fp_rdm_r = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_fp_rdm_pooled_r.wig'
mrna_rdm_f = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_mrna-rdm-pooled_f.wig'
mrna_rdm_r = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_mrna-rdm-pooled_r.wig'
fp_min_f = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Unused/GSE53767_fp_min_pooled_f.wig'
fp_min_r = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Unused/GSE53767_fp_min_pooled_r.wig'
geneAnnotationDir = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/CDS_MG1655.mochiview.txt'


def getStartCodons():
    startTime = time.time()
    startCodons = {}

    #Get the dictionary of genes from the mochi file.
    fileName = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\RNA Seq\Jing\WIGFiles\Used\CDS_MG1655.mochiview.txt"
    geneDict = Reads.makeDictFromMochi(fileName)

    #Turn the genome into a string.
    genomelocation = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/MG1655NC_000913.fasta'
    fullseq = Reads.convertGenomeToString(genomelocation)
    print "finished loading genome"
    print len(fullseq), "nucleotides"

    for key in geneDict:
        gene = geneDict[key]

        if gene["STRAND"]=="+":
            #continue
            start = gene["START"]
            #adjust because indexing of fullseq starts at 0.
            startCodon = fullseq[start-1:start+2]
        elif gene["STRAND"]=="-":
            #continue
            start = gene["END"]
            #reverse complement
            startCodon = fullseq[start-1] + fullseq[start-2] + fullseq[start-3]
            startCodon = startCodon.replace("A","x")
            startCodon = startCodon.replace("T","A")
            startCodon = startCodon.replace("x","T")
            startCodon = startCodon.replace("G","x")
            startCodon = startCodon.replace("C","G")
            startCodon = startCodon.replace("x","C")
##        if (gene["END"]-gene["START"]+1)%3 !=0:
##            print gene["GENE_NAME"]
##            print key, gene
##            print startCodon
##            pause = raw_input("pause")
        if startCodon not in startCodons:
            startCodons[startCodon] = 1
        else:
            startCodons[startCodon] = startCodons[startCodon] + 1
    endTime = time.time()
    print "This took",(endTime-startTime),"seconds."
    return startCodons

def getCodonFrequency(genomelocation, files, strands, pairings, geneAnnotationDir):
    startTime = time.time()
    
    allGeneDicts = Reads.wrapper(genomelocation, files, strands, pairings, geneAnnotationDir)
    allCodonUsages = []
    while len(allCodonUsages)<len(allGeneDicts):
        allCodonUsages.append({})
    #Make the usage dictionary for the ribo-seq data
    for i in range(len(allCodonUsages)):
        codonUsage = allCodonUsages[i]
        geneDict = allGeneDicts[i]
        for key in geneDict:
            gene = geneDict[key]
            if gene["GENE_NAME"]=="prfB" or gene["GENE_NAME"]=="dnaX" or gene["GENE_NAME"]=="FdhF" or gene["GENE_NAME"]=="FdoG"\
               or gene["GENE_NAME"]=="FdnG" or gene["GENE_NAME"]=="ArfA" or gene["GENE_NAME"]=="SecM" or gene["GENE_NAME"]=="TnaC":
                continue
            seq = gene["SEQ"]
            scores = gene["SEPARATE_SCORE"]
            for i in range(len(seq)/3):
                codon = seq[i*3:i*3+3]
                codonScore = scores[i*3]+scores[i*3+1]+scores[i*3+2]
                if codon in codonUsage:
                    codonUsage[codon][0] += 1
                    codonUsage[codon][1] += codonScore
                else:
                    codonUsage[codon] = [1,codonScore]
        print 100.0*i/len(allCodonUsages),"percent complete with making codon usage dictionaries."
    endTime = time.time()
    print "This full process took",(endTime-startTime),"seconds."
    return allCodonUsages, allGeneDicts

def makePlotForRDMCodonUsage(genomelocation='', files=[], strands=[], pairings=[], geneAnnotationDir=''):
    #Values specific to this test
    genomelocation = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/MG1655NC_000913.fasta'
    files = [fp_rdm_f,fp_rdm_r,mrna_rdm_f,mrna_rdm_r]
    strands = ['+','-','+','-']
    pairings = [(0,1),(2,3)]
    geneAnnotationDir = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/CDS_MG1655.mochiview.txt'
    
    allCodonUsages, allGeneDicts = getCodonFrequency(genomelocation, files, strands, pairings, geneAnnotationDir)
    codonUsage1 = allCodonUsages[0]
    codonUsage2 =allCodonUsages[1]
    keys = codonUsage1.keys()
    keys.sort()
    combined = {}

    for codon in keys:
        if codonUsage1[codon][0] == codonUsage2[codon][0]:
            combined[codon] = [codonUsage1[codon][0],codonUsage1[codon][1]/codonUsage2[codon][1]]
        else:
            print codon, "There is something wrong with this codon."
    
    maxCodonUsage = 0
    sumOfScore = 0

    for codon in keys:
        codonInfo = combined[codon]
        sumOfScore += codonInfo[1]
        if codonInfo[0]>maxCodonUsage:
            maxCodonUsage = codonInfo[0]
    avgScore = sumOfScore/len(combined)
    print avgScore, "is the average score."
    scaledCombined = {}
    for codon in keys:
        scaledCombined[codon] = [combined[codon][0]/float(maxCodonUsage),combined[codon][1]/avgScore]
    
    xVals = []
    yVals = []
    for codon in keys:
        xVals.append(scaledCombined[codon][0])
        yVals.append(scaledCombined[codon][1])
        
    plt.figure()
    plt.xlim(0,1)
    plt.ylim(0,3)
    plt.plot(xVals, yVals, 'ro')
    plt.xlabel('Codon Usage (AU)')
    plt.ylabel('Ribosome Occupancy (AU)')
    plt.title('Codons')

def compareMinToRich(genomelocation = "", files = [], strands = [], pairings = [], geneAnnotationDir = ""):
    
    #Comment this out when you want to specify which to use
    #Values specific to this test
    genomelocation = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/MG1655NC_000913.fasta'
    files = [fp_rdm_f,fp_rdm_r,fp_min_f,fp_min_r]
    strands = ['+','-','+','-']
    pairings = [(0,1),(2,3)]
    geneAnnotationDir = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/CDS_MG1655.mochiview.txt'
    
    allCodonUsages, allGeneDicts = getCodonFrequency(genomelocation, files, strands, pairings, geneAnnotationDir)
    codonUsage1 = allCodonUsages[0]
    codonUsage2 = allCodonUsages[1]
    geneDict1 = allGeneDicts[0]
    geneDict2 = allGeneDicts[1]
    codonXVals = []
    codonYVals = []
    geneXVals = []
    geneYVals = []

    codonKeys = codonUsage1.keys()
    geneKeys = geneDict1.keys()

    for codon in codonKeys:
        codonYVals.append(codonUsage1[codon][1])
        codonXVals.append(codonUsage2[codon][1])
    for geneID in geneKeys:
        geneYVals.append(geneDict1[geneID]['SCORE'])
        geneXVals.append(geneDict2[geneID]['SCORE'])
    
    plt.figure()
    plt.plot(codonXVals, codonYVals, 'ro')
    plt.xlabel('Minimal Media')
    plt.ylabel('Rich Defined Media')
    plt.title('Codon Score Comparison')

    plt.figure()
    plt.plot(geneXVals, geneYVals, 'ro')
    plt.xlabel('Minimal Media')
    plt.ylabel('Rich Defined Media')
    plt.title('Gene Score Comparison')

