##These functions give the truth value and the gate count.

##This needs to be generalized to work for varying inputs.
##This may also need to be updated when a representation for NAND, NOR, IMPLIES
##and NIMPLIES is decided on.
def getTruthValue(circuit, numInputs):
    """
        Takes in a circuit and number of inputs and returns the string
        representation of the truth value the circuit evaluates to
    """
    if numInputs == 3:
        circuit = circuit.replace("0","00000000")
        circuit = circuit.replace("1","11111111")
        circuit = circuit.replace("a","00001111")
        circuit = circuit.replace("b","00110011")
        circuit = circuit.replace("c","01010101")
    if numInputs == 2:
        circuit = circuit.replace("0","0000")
        circuit = circuit.replace("1","1111")
        circuit = circuit.replace("a","0011")
        circuit = circuit.replace("b","0101")
    indeces = findSmallestParentheses(circuit)
    while indeces != None:
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
        #In case we run into (a) or ((a.b)) or ()
        else:
            answer = subCircuit[1:-1]
        circuit = circuit.replace(subCircuit,answer)
        indeces = findSmallestParentheses(circuit)
    return circuit

##This might be updated when I decide on the appropriate representation for
##NAND, NOR, IMPLIES, and NIMPLIES.
def circuitCost(circuit, cost={"~":1,"&":1,"@":1,"+":1,"^":1,".":1,"=":1,">":1,"$":1}):
    """
        Takes in a circuit and a dictionary that maps each logic gate to it's
        cost and returns the cost of a circuit.
        Defaults to each logic operator costing 1.
    """
    operators = cost.keys()
    operatorCount = {}
    for operator in operators:
        operatorCount[operator] = 0
    indeces = findSmallestParentheses(circuit)
    if indeces == None:
        return 0
    count = 1
    while indeces != None:
##        print circuit
        subCircuit = circuit[indeces[0]:indeces[1]+1]
        for operator in operators:
            operatorIndex = subCircuit.find(operator)
            if operatorIndex != -1:
                operatorCount[operator] = operatorCount[operator] + 1
                #This will make sure we don't count (a+b) and (b+a) separately

                if operator == "&" or operator == "@" or operator == "+" \
                   or operator == "^" or operator == "." or operator == "=":
                    subCircuit2 = "("+subCircuit[operatorIndex+1:-1]+operator+subCircuit[1:operatorIndex]+")"
                    circuit = circuit.replace(subCircuit2,str(count))
                break
        
        circuit = circuit.replace(subCircuit,str(count))
        indeces = findSmallestParentheses(circuit)
        count += 1
    
##    print circuit
    totalCost = 0
    for operator in operators:
        totalCost += cost[operator]*operatorCount[operator]
    return totalCost


##Logic Operator Functions

def NOT(s):
    """
        Takes a truth value enclosed in parentheses with a NOT symbol returns
        the truth value after performing the NOT
    """
    operatorIndex = s.index("~")
    #get the two truth values removing the parentheses
    a = s[operatorIndex+1:-1]
    answer = ""
    for i in xrange(len(a)):
        if a[i]=="0":
            answer += "1"
        else:
            answer += "0"
    return answer

def AND(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        an AND symbol between them and returns the truth value after performing
        the AND
    """
    operatorIndex = s.index("&")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]=="0" or b[i]=="0"):
            answer += "0"
        else:
            answer += "1"
    return answer

##The symbol for this might change
def NAND(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        a NAND symbol between them and returns the truth value after performing
        the NAND
    """
    operatorIndex = s.index("@")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]=="0" or b[i]=="0"):
            answer += "1"
        else:
            answer += "0"
    return answer

def OR(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        an OR symbol between them and returns the truth value after performing
        the OR
    """
    operatorIndex = s.index("+")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]=="1" or b[i]=="1"):
            answer += "1"
        else:
            answer += "0"
    return answer

def XOR(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        an XOR symbol between them and returns the truth value after performing
        the XOR
    """
    operatorIndex = s.index("^")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]==b[i]):
            answer += "0"
        else:
            answer += "1"
    return answer

##The symbol for this might change
def NOR(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        a NOR symbol between them and returns the truth value after performing
        the NOR
    """
    operatorIndex = s.index(".")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]=="1" or b[i]=="1"):
            answer += "0"
        else:
            answer += "1"
    return answer

def XNOR(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        an XNOR symbol between them and returns the truth value after performing
        the XNOR
    """
    operatorIndex = s.index("=")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]==b[i]):
            answer += "1"
        else:
            answer += "0"
    return answer

##The symbol for this might change
def IMPLIES(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        an IMPLIES symbol between them and returns the truth value after
        performing the IMPLIES
    """
    operatorIndex = s.index(">")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]=="0" or b[i]=="1"):
            answer += "1"
        else:
            answer += "0"
    return answer

##The symbol for this might change
def NIMPLIES(s):
    """
        Takes two truth values of the same length enclosed in parentheses with
        a NIMPLIES symbol between them and returns the truth value after
        performing the NIMPLIES
    """
    operatorIndex = s.index("$")
    #get the two truth values removing the parentheses
    a = s[1:operatorIndex]
    b = s[operatorIndex+1:-1]
    if len(a)!=len(b):
        print "The values are different sizes. Can't perform operation."
        return None
    answer = ""
    for i in xrange(len(a)):
        if (a[i]=="0" or b[i]=="1"):
            answer += "0"
        else:
            answer += "1"
    return answer

##Buffer incase needed
def BUFFER(s):
    """
        Takes a truth value and returns the truth value
    """
    return s


## Helper Functions
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
        print truth,":",foundCircuits[truth]

##Needs to be generalized for circuits that are not three inputs
def containsEqualSubCircuit(circuit, numInputs):
    """
        Checks if the circuit has a subcircuit in it with the same truth value
    """
    truthValue = getTruthValue(circuit,numInputs)
    if numInputs == 3:
        circuit = circuit.replace("0","00000000")
        circuit = circuit.replace("1","11111111")
        circuit = circuit.replace("a","00001111")
        circuit = circuit.replace("b","00110011")
        circuit = circuit.replace("c","01010101")
    if numInputs == 2:
        circuit = circuit.replace("0","0000")
        circuit = circuit.replace("1","1111")
        circuit = circuit.replace("a","0011")
        circuit = circuit.replace("b","0101")
    indeces = findSmallestParentheses(circuit)
    while indeces != None:
        if circuit.find(truthValue) != -1:
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
