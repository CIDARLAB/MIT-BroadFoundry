##These are the functions that are used for saving dictionaries of truth values and
##their minimal circuits to text files and retrieving them.

def findLeftSingleQuotes(s):
    """
        Takes in a string that contains a single quote.
        Returns the indeces of the first pair of quotes.
        Returns None if there are no quotes.
    """
    try:
        startIndex = s.index("'")
        endIndex = s.index("'",startIndex+1)
        return (startIndex,endIndex)
    except ValueError:
        return None

def writeToFile(foundCircuits, fileName):
    """
        Takes a dictionary of truth values mapped to their minimum circuits and
        writes it to a file with the given fileName
    """
    myFile = open(fileName,'w')
    keys = foundCircuits.keys()
    keys.sort()
    for key in keys:
        s = key + " " + str(foundCircuits[key]) + "\n"
        myFile.write(s)
    myFile.close()

def getFromFile(fileName):
    """
        Opens a file and creates a dictionary of truth values mapped to a list of
        their minimum gate circuits
    """
    foundCircuits = {}
    myFile = open(fileName,'r')
    for line in myFile:
        index = line.index(" ")
        key = line[0:index]
        unformattedValue = line[index:]
        value = []
        
        #We are going to take all the circuits from the unformatted value and
        #append it to value
        quoteIndeces = findLeftSingleQuotes(unformattedValue)
        while quoteIndeces != None:
            #Add the circuit to value and remove it from the string
            value.append(unformattedValue[quoteIndeces[0]+1:quoteIndeces[1]])
            unformattedValue = unformattedValue[:quoteIndeces[0]] + unformattedValue[quoteIndeces[1]+1:]
            quoteIndeces = findLeftSingleQuotes(unformattedValue)
        foundCircuits[key] = value
    myFile.close()
    print str(len(foundCircuits)),"truth values found"
    return foundCircuits
