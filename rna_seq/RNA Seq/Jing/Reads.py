import time
import copy
import matplotlib.pyplot as plt
import numpy as np
import scipy

#File Locations
genomelocation = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/MG1655NC_000913.fasta'
fp_rdm_f = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_fp_rdm_pooled_f.wig'
fp_rdm_r = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_fp_rdm_pooled_r.wig'
mrna_rdm_f = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_mrna-rdm-pooled_f.wig'
mrna_rdm_r = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/GSE53767_mrna-rdm-pooled_r.wig'
fp_min_f = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Unused/GSE53767_fp_min_pooled_f.wig'
fp_min_r = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Unused/GSE53767_fp_min_pooled_r.wig'
geneAnnotationDir = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/CDS_MG1655.mochiview.txt'

def getReads(location, fullseq):
    """
        Takes the file location for a wig file and a genome as a string.
        Returns a dictionary of genome spots mapped to a list containing the letter at that spot and the scaled reads.
    """
    #Open the file, initialize the dictionary, and get the length of the genome.
    f = open(location, 'r')
    reads = {}
    genomeSize = len(fullseq)

    #Go through the wig file and add to the dictionary the positions in the file and the nuceotide at that position and the number of reads.
    for line in f:
        parts = line.split("\t")
        reads[int(parts[0])] = [fullseq[int(parts[0])-1],float(parts[1])]
    f.close()
    #For anything that wasn't in the file, add it to the dictionary with its nucleotide and 0 reads.
    for i in range(1,genomeSize+1):
        if i not in reads:
            reads[i] = [fullseq[i-1],0.0];
    return reads

def convertGenomeToString(genomelocation):
    genome = open(genomelocation, 'r')
    fullseq = ""
    genome.readline()
    for line in genome:
        fullseq += line[:-1]
    genome.close()
    return fullseq

def getAllReads(genomelocation, files, strands, pairings):
    """
        Takes a string that is the directory of the genome, a list of strings that are the
        directories for each of the wig files, and a list of strings that define each file to be forward or reverse
        using + or - in the same order as the directories in files. Also takes in pairings which is a list of tuples
        that indicates the indeces paired together for reverse forward pairs.
        Example:
        files = [fp_rdm_f,fp_min_f,fp_rmd_r,fp_min_r]
        starnds = ['+','+','-','-']
        pairings = [(0,2),(1,3)]
        Returns a list of dictionaries produced using getReads, the list strands, the list of tuples pairings, and the genome as a string.
    """

    #Start the timer
    #startTime = time.time()

    #Turn the genome into a string.
    print "loading genome"
    fullseq = convertGenomeToString(genomelocation)
    print "finished loading genome"
    print len(fullseq), "nucleotides"
    
    #Use the string form of the genome along with the wig files to produce each of the dictionaries.
    print "Making Dictionaries of reads"
    reads = []
    count = 0
    for fileLoc in files:
        reads.append(getReads(fileLoc,fullseq))
        count += 1
        print 100*float(count)/len(files),"% complete"

##    #Stop the timer and print time.
##    endTime = time.time()
##    print "This took", endTime-startTime,"seconds."
    
    return reads, strands, pairings, fullseq

def makeDictFromMochi(fileName):
    """
        Takes the file location of the gene anotation file and splits it into a dictionary where the key is the b#### feature name, and the value is
        a dictionary of all the remaining properties of that gene.
    """

    #Open the file.
    myFile = open(fileName, 'r')
    
    #The first line is jut text that will be used as keys for dictionary in the values.
    keys = myFile.readline().split('\t')

    #Remove whitespace.
    for i in range(len(keys)):
        keys[i] = keys[i].strip()
    
    #Find the feature name index and initialize the empty dictionary.
    geneIDIndex = keys.index("FEATURE_NAME")
    geneDict = {}

    for line in myFile:
        #For the remaining lines with the actual content split the values up and take out the feature name and use it as a key and initialize the dictionary.
        geneContent = line.split('\t')
        geneID = geneContent[geneIDIndex]
        geneDict[geneID] = {}
        #For all other entries put them in the initialied dictionary for that gene with the tag from the first line of the document as their key.
        for i in range(len(keys)):
            if i!=geneIDIndex:
                try:
                    geneDict[geneID][keys[i]] = int(geneContent[i].strip())
                except ValueError:
                    geneDict[geneID][keys[i]] = geneContent[i].strip()
    myFile.close()                    
    return geneDict

def updateDictWithCountAndSeq(geneDict,reads,sign):
    """
        This takes in a dictionary produced by makeDictFromMochi made from the gene annotation file, and a dictionary made from findReads() made from the wig file.
        It also takes in whether the wig file was the forward or reverse sequence. If forward, sign=="+" if reverse sign =="-".
        Updates geneDict to include the sequence of the coding region of the gene, the total read count, and an ordered list containing
        the read  count for each nucleotide.
    """

    #get the list of genes and determine whether or not we will need to do the reverse complement.
    genelist = geneDict.keys()
    if sign == "-":
        complement = True
    else:
        complement = False

    #go through each gene and get the gene information.
    for geneName in genelist:
        gene = geneDict[geneName]
        #only continue if the gene strand matches the sign.
        if gene['STRAND'] == sign:
            #get the starting and ending nucleotide and initialize the values for the reads and sequence
            geneStart = gene['START']
            geneEnd = gene['END']
            seq = ""
            score = 0
            scoreList = []

            #get the ragnge of nucleotide numbers we will go through.
            nucNum = range(geneStart,geneEnd+1)

            #reverse the gene if necessary
            if complement:
                    nucNum.reverse()
                    
            for i in nucNum:
                #complement the gene if necessary
                if complement:
                    nuc = reads[i][0]
                    if nuc=="A":
                        seq+="T"
                    elif nuc=="T":
                        seq+="A"
                    elif nuc=="C":
                        seq+="G"
                    elif nuc=="G":
                        seq+="C"
                else:
                    seq += reads[i][0]
                score += reads[i][1]
                scoreList.append(reads[i][1])

            #add it to the properties of the gene. note that this mutates (changes) geneDict so nothing needs to be returned.
            gene['SEQ'] = seq
            gene['SCORE'] = score
            gene['SEPARATE_SCORE'] = scoreList

def wrapper(genomelocation, files, strands, pairings, geneAnnotationDir):
    """
        Wraps together all of the above functions into one. Returns a list of dictionaries with the genes b#### feature name as the keys and a dictionary of all of its properties,
        including the sequence of the translated region and the read counts, as the values. The first dictionary uses the fp wig files to get the read counts, the
        second uses the mrna wig files to get the read counts.
    """
    if len(files)%2!=0 or strands.count('+')!=strands.count('-') or len(pairings)!=len(files)/2:
        print "Make sure you have included only one forward and reverse file for each pairing and that pairings were done appropriately"
        return
    
    #startTime = time.time()

    #get the dictionaries from the wig files.
    reads, strands, pairings, fullseq = getAllReads(genomelocation, files, strands, pairings)

    #Make the dictionaries from the gene annotation file. Make a deepcopy of it for each pair of wig files.
    allGeneDicts = [makeDictFromMochi(geneAnnotationDir)]
    while len(allGeneDicts)<len(pairings):
        allGeneDicts.append(copy.deepcopy(allGeneDicts[0]))
    print "Finished making template dicts"

    #Update the geneDict for all the + genes then for all the - genes.
    for i in range(len(allGeneDicts)):
        geneDict = allGeneDicts[i]
        currPairing = pairings[i]
        currPairing[0]
        currPairing[1]
        if strands[currPairing[0]] == strands[currPairing[1]]:
            print "Something is wrong with the strands list."
            return
        updateDictWithCountAndSeq(geneDict,reads[currPairing[0]],strands[currPairing[0]])
        print 100.0*(i+0.5)/len(allGeneDicts),"% complete updating dictionaries with gene sequences and counts."
        updateDictWithCountAndSeq(geneDict,reads[currPairing[1]],strands[currPairing[1]])
        print 100.0*(i+1.0)/len(allGeneDicts),"% complete updating dictionaries with gene sequences and counts."
    print "Done updating all dicts."
    
    return allGeneDicts

def remakewig(wigFile, genomeSize=4639675):
    """
        Takes in a wig file and remakes the file so that it is in a format that can be used in MochiView. Genome size of E.Coli==4639675
    """
    startTime= time.time()
    #This is how a file must start to be recognizable by MochiView
    s = "track type=wiggle_0\nvariableStep chrom=NC_000913_2\n"

    #Initialize the dictionary. Open the old wig file and create and open the new wig file.
    allPoints = {}
    myFile = open(wigFile, 'r')
    myFile2 = open(wigFile[:-4]+"_editted.wig",'w')
    
    print "getting existing nucleotides"
    for line in myFile:
        #Split the line and take the first value to be the nucleotide number and the second to be the read count.
        try:
            x = line.split('\t')
            allPoints[int(x[0])] = x[1].strip()
        #in case the line is not in the proper format
        except ValueError:
            pass

    #Write it to the file
    print "starting to write"
    for i in range(1,genomeSize+1):
        #If the nucleotide number was in the old file, use it. Otherwise put in a value of 0.0
        try:
            s = s+str(i)+'\t'+str(allPoints[i])+'\n'
        except KeyError:
            s = s+str(i)+'\t'+'0.0'+'\n'
        #Every time s gets to this size, write it to a file and reset s.
        if len(s)>=500000:
            #Print progress
            print round((100.0*i/4639675.0),2),"% completed"
            myFile2.write(s)
            s = ""
    #Write whatever is left over to the file.
    myFile2.write(s)        
    print "done"

    #Close the files.
    myFile.close()
    myFile2.close()
    endTime = time.time()
    print "This took",(endTime-startTime),"seconds."
    return

def compareFirstSecondHalf():
    files = [fp_rdm_f,fp_rdm_r]

    strands = ['+','-']

    pairings = [(0,1)]

    x = wrapper(genomelocation,files,strands,pairings,geneAnnotationDir)
    xValsRed = []
    yValsRed = []
    xValsBlack = []
    yValsBlack = []
    skip = ['prfB','dnaX','FdhF','FdoG','FdnG','SecM','TnaC']
    for key in x[0]:
        if x[0][key]['GENE_NAME'] in skip:
            continue
        elif len(x[0][key]['SEQ'])<60 or sum(x[0][key]['SEPARATE_SCORE'][15:len(x[0][key]['SEPARATE_SCORE'])/2])<64 or sum(x[0][key]['SEPARATE_SCORE'][len(x[0][key]['SEPARATE_SCORE'])/2:-15])<64:
            firstHalf = x[0][key]['SEPARATE_SCORE'][15:len(x[0][key]['SEPARATE_SCORE'])/2]
            secondHalf = x[0][key]['SEPARATE_SCORE'][len(x[0][key]['SEPARATE_SCORE'])/2:-15]
            xValsBlack.append(sum(firstHalf)/len(firstHalf))
            yValsBlack.append(sum(secondHalf)/len(secondHalf))
        else:
            firstHalf = x[0][key]['SEPARATE_SCORE'][15:len(x[0][key]['SEPARATE_SCORE'])/2]
            secondHalf = x[0][key]['SEPARATE_SCORE'][len(x[0][key]['SEPARATE_SCORE'])/2:-15]
            xValsRed.append(sum(firstHalf)/len(firstHalf))
            yValsRed.append(sum(secondHalf)/len(secondHalf))
    
    
    coeff = rsquared(xValsRed,yValsRed)
    predicted = []
    for i in range(len(xValsRed)):
        predicted.append(coeff[0]*xValsRed[i])# + coeff[1])
    print "line equation:", ("y="+str(coeff[0])+"x")#+"+str(coeff[1]))
    print "r^2: ", coeff[2]
    
    plt.figure()
    plt.plot(xValsBlack, yValsBlack, 'k.', ms = 2.0)
    plt.plot(xValsRed, yValsRed, 'r.',xValsRed,predicted, 'b.',ms=2.0)
    plt.xlim(10**-3,10**5)
    plt.ylim(10**-3,10**5)
    plt.xscale('log')
    plt.yscale('log')    
    
    plt.xlabel('Ribosome density of 2nd half')
    plt.ylabel('Ribosome density of 1st half')
##    plt.title('Codons')
    


def rsquared(x,y):
    ''' Returns the r^2 value '''
    slope,intercept,r_value,p_value,std_err = scipy.stats.linregress(x,y)
    return slope, intercept, r_value**2

def getStartCodons():
    startTime = time.time()
    startCodons = {}

    #Get the dictionary of genes from the mochi file.
    fileName = "C:\Users\Arinze\SkyDrive\UROP_Summer_2015\RNA Seq\Jing\WIGFiles\Used\CDS_MG1655.mochiview.txt"
    geneDict = makeDictFromMochi(fileName)

    #Turn the genome into a string.
    genomelocation = 'C:/Users/Arinze/SkyDrive/UROP_Summer_2015/RNA Seq/Jing/WIGFiles/Used/MG1655NC_000913.fasta'
    fullseq = convertGenomeToString(genomelocation)
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
    
    allGeneDicts = wrapper(genomelocation, files, strands, pairings, geneAnnotationDir)
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


