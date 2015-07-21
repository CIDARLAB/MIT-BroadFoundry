from Truths_and_Gates import *
##These functions are functions useful for finding truth values and gate counts

def toBinary(num,numInputs):
    """
        Takes in a decimal number and the number of inputs for the circuit
        and returns a binary representation as a string for it of the appropriate length
        given the number of inputs.
    """
    L = 1
    for i in range(numInputs):
        L*=2
    s = "{0:b}".format(num)
    while len(s)<L:
        s = "0"+s
    return s

def findLeftParentheses(circuit):
    """
        Finds the indeces of the leftmost open parenthesis and its matching
        close parenthesis returns their indeces
    """
    try:
        startIndex = circuit.index("(")
        L = len(circuit)
        counter = 1
        for i in xrange(startIndex+1,L):
            if circuit[i] == "(":
                counter += 1
            elif circuit[i] == ")":
                counter -= 1
                if counter == 0:
                    endIndex = i
                    break
        return (startIndex, endIndex)
    except ValueError:
        return None
    except NameError:
        print "these parentheses aren't properly matched"
        return None

def findSmallestParentheses(circuit):
    """
        Returns the indeces of a pair of parentheses that do not contain
        any parentheses in it
    """
    L = len(circuit)
    for i in xrange(L):
        if circuit[i]=="(":
            startIndex = i
        if circuit[i]==")":
            endIndex = i
            try:
                return (startIndex,endIndex)
            except NameError:
                print "parentheses are not properly matched"
                return None
    return None

def sortAndPrint(foundCircuits):
    """
        Takes a ditionary of truth values mapped to a list of their circuits
        and prints each truth value and the list of its minimum circuits in
        order of truth values
    """
    truths = foundCircuits.keys()
    truths.sort()
    for truth in truths:
        print truth," : ",foundCircuits[truth]

#Needs to be generalized for circuits that are not three inputs
def containsEqualSubCircuit(circuit, numInputs):
    """
        Checks if the circuit has a subcircuit in it with the same truth value
    """
    truthValue = getTruthValue(circuit,numInputs)
    circuit = circuit.replace("0","00000000")
    circuit = circuit.replace("1","11111111")
    circuit = circuit.replace("a","00001111")
    circuit = circuit.replace("b","00110011")
    circuit = circuit.replace("c","01010101")
    indeces = findSmallestParentheses(circuit)
    while indeces != None:
        #If I used != -1 this would print true for circuit that are (a)
        #which we don't want even if we don't expect this ciruit to come up
        if circuit.find(truthValue) > 0:
            return True
        subCircuit = circuit[indeces[0]:indeces[1]+1]
        if subCircuit.find("~") != -1:
            answer = NOT(subCircuit)
        elif subCircuit.find("&") != -1:
            answer = AND(subCircuit)
        elif subCircuit.find("@") != -1:
            answer = NAND(subCircuit)
        elif subCircuit.find("+") != -1:
            answer = OR(subCircuit)
        elif subCircuit.find("^") != -1:
            answer = XOR(subCircuit)
        elif subCircuit.find(".") != -1:
            answer = NOR(subCircuit)
        elif subCircuit.find("=") != -1:
            answer = XNOR(subCircuit)
        elif subCircuit.find(">") != -1:
            answer = IMPLIES(subCircuit)
        elif subCircuit.find("$") != -1:
            answer = NIMPLIES(subCircuit)
        else:
            answer = subCircuit[1:-1]
        circuit = circuit.replace(subCircuit,answer)
        indeces = findSmallestParentheses(circuit)
    return False

