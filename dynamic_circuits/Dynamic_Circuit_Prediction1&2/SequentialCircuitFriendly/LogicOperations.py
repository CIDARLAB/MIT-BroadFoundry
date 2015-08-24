def NOT(tvList):
    """
    Takes in a list of length one that contains a single value 1 or 0
    as a string returns the opposite of the value
    """
    assert len(tvList)==1 and type(tvList)==list

    if tvList[0] == "1":
        return "0"
    return "1"
    
def BUF(tvList):
    """
    Takes in a list of length one that contains a single value 1 or 0
    as a string returns the opposite of the value
    """
    assert len(tvList)==1 and type(tvList)==list

    if tvList[0] == "0":
        return "0"
    return "1"
    
def NOR(tvList):
    assert len(tvList)>=2 and type(tvList)==list
    for i in tvList:
        if i=="1":
            return "0"
    return "1"
    
def OR(tvList):
    assert len(tvList)>=2 and type(tvList)==list
    for i in tvList:
        if i=="1":
            return "1"
    return "0"
    
def AND(tvList):
    assert len(tvList)>=2 and type(tvList)==list
    for i in tvList:
        if i=="0":
            return "0"
    return "1"
    
def NAND(tvList):
    assert len(tvList)>=2 and type(tvList)==list
    for i in tvList:
        if i=="0":
            return "1"
    return "0"
    
def XOR(tvList):
    assert len(tvList)>=2 and type(tvList)==list
    state = False
    for i in tvList:
        if i=="1" and state==False:
            state=True
        #if we see a second true then return 0
        elif i=="1" and state==True:
            return "0"
    if state:
        return "1"
    return "0"
    
def XNOR(tvList):
    """
    Same as XOR except returns 0 when XOR would return 1 and 
    returns 1 when XOR would return 0.
    """
    assert len(tvList)>=2 and type(tvList)==list
    state = False
    for i in tvList:
        if i=="1" and state==False:
            state=True
        elif i=="1" and state==True:
            return "1"
    if state:
        return "0"
    return "1"