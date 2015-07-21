import inputs

#Change Functions
def getmRNAChange(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    #Determine what type of gate it is and pass it to the appropriate function.
    if gateProperties['TYPE']=='NOR':
        return NOR(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein)
    elif gateProperties['TYPE']=='NOT':
        return NOT(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein)
    elif gateProperties['TYPE']=='OR':
        return OR(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein)
    elif gateProperties['TYPE']=='BUFFER':
        return BUFFER(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein)
    #In case the gate type is not recognized.
    else:
        print "This is not a recognized gate type."
        return


def getProteinChange(initProteinVal,initmRNAVal,gateProperties):
    #dX_dt = pB*mX - ap*Xi
    return gateProperties['pB']*initmRNAVal - gateProperties['ap']*initProteinVal



#Gates
def NOR(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    '''
        Given appropriate inputs will do a double repression and calculate the
        change in mRNA value.
    '''
    #dmX_dt = -am*mX + mBx1*(Km1**n1)/(Input1**n1+Km1**n1) + mB2*(Km2**n2)/(Input2**n2+Km2**n2)
    #First term
    firstTerm = -gateProperties['am']*initmRNAVal
    
    #Second term.
    #Getting values from the first input
    mB1 = gateProperties['mB1']
    Km1 = gateProperties['Km1']
    n1 = gateProperties['n1']
    InputName1 = gateProperties['INPUT1']
    InputProperties1 = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName1)]
    #If this is another gate then just get the concentration of the protein from it.
    if InputProperties1['TYPE'] != 'INPUT':
        Input1 = initProtein[logic_gate_names.index(InputName1)]
    #If it is an input, then use the equation specified for the input to determine
    #the concentration of the input at that time step.
    else:
        Inputfunc1 = InputProperties1['INPUT']
        Input1 = Inputfunc1[0](t,*Inputfunc1[1:])
    secondTerm = repressor(mB1,Km1,n1,Input1)
    
    #Third Term
    #Getting values from the second input
    mB2 = gateProperties['mB2']
    Km2 = gateProperties['Km2']
    n2 = gateProperties['n2']
    InputName2 = gateProperties['INPUT2']
    InputProperties2 = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName2)]
    #If this is another gate then just get the concentration of the protein from it.
    if InputProperties2['TYPE'] != 'INPUT':
        Input2 = initProtein[logic_gate_names.index(InputName2)]
    #If it is an input, then use the equation specified for the input to determine
    #the concentration of the input at that time step.
    else:
        Inputfunc2 = InputProperties2['INPUT']
        Input2 = Inputfunc2[0](t,*Inputfunc2[1:])
    thirdTerm = repressor(mB2,Km2,n2,Input2)
    return firstTerm + secondTerm + thirdTerm


def NOT(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    '''
        Given the necessary information will perform a simple repression by the input
    '''
    #dm_dt = -am*mX + mB*(Km**n)/(Input**n+Km**n)
    #First term
    firstTerm = -gateProperties['am']*initmRNAVal

    #Second term
    mB = gateProperties['mB']
    Km = gateProperties['Km']
    n = gateProperties['n']
    InputName = gateProperties['INPUT']
    InputProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName)]
    #If this is another gate then just get the concentration of the protein from it.
    if InputProperties['TYPE'] != 'INPUT':
        Input = initProtein[logic_gate_names.index(InputName)]
    #If it is an input, then use the equation specified for the input to determine
    #the concentration of the input at that time step.
    else:
        Inputfunc = InputProperties['INPUT']
        Input = Inputfunc[0](t,*Inputfunc[1:])
    secondTerm = repressor(mB,Km,n,Input)
    
    return firstTerm + secondTerm

    
def OR(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    '''
        Given appropriate inputs will do a double activator and calculate the
        change in mRNA value.
    '''
    #dmX_dt = -am*mX + mB1*(Input1**n1)/(Input1**n1+Km1**n1) + mB2*(Input2**n2)/(Input2**n2+Km2**n2)
    #First term
    firstTerm = -gateProperties['am']*initmRNAVal
    
    #Second term
    mB1 = gateProperties['mB1']
    Km1 = gateProperties['Km1']
    n1 = gateProperties['n1']
    InputName1 = gateProperties['INPUT1']
    InputProperties1 = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName1)]
    #If this is another gate then just get the concentration of the protein from it.
    if InputProperties1['TYPE'] != 'INPUT':
        Input1 = initProtein[logic_gate_names.index(InputName1)]
    #If it is an input, then use the equation specified for the input to determine
    #the concentration of the input at that time step.
    else:
        Inputfunc1 = InputProperties1['INPUT']
        Input1 = Inputfunc1[0](t,*Inputfunc1[1:])
    secondTerm = activator(mB1,Km1,n1,Input1)
    
    #Third Term
    mB2 = gateProperties['mB2']
    Km2 = gateProperties['Km2']
    n2 = gateProperties['n2']
    InputName2 = gateProperties['INPUT2']
    InputProperties2 = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName2)]
    #If this is another gate then just get the concentration of the protein from it.
    if InputProperties2['TYPE'] != 'INPUT':
        Input2 = initProtein[logic_gate_names.index(InputName2)]
    #If it is an input, then use the equation specified for the input to determine
    #the concentration of the input at that time step.
    else:
        Inputfunc2 = InputProperties2['INPUT']
        Input2 = Inputfunc2[0](t,*Inputfunc2[1:])
    thirdTerm = activator(mB2,Km2,n2,Input2)
    return firstTerm + secondTerm + thirdTerm

def BUFFER(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    '''
        Given the necessary information will perform a simple activator by the input
    '''
    #dm_dt = -am*mX + mB*(Input**n)/(Input**n+Km**n)
    #First term
    firstTerm = -gateProperties['am']*initmRNAVal

    #Second term
    mB = gateProperties['mB']
    Km = gateProperties['Km']
    n = gateProperties['n']
    InputName = gateProperties['INPUT']
    InputProperties = input_and_logic_gate_dictionaries[input_and_logic_gate_names.index(InputName)]
    #If this is another gate then just get the concentration of the protein from it.
    if InputProperties['TYPE'] != 'INPUT':
        Input = initProtein[logic_gate_names.index(InputName)]
    #If it is an input, then use the equation specified for the input to determine
    #the concentration of the input at that time step.
    else:
        Inputfunc = InputProperties['INPUT']
        Input = Inputfunc[0](t,*Inputfunc[1:])
    secondTerm = activator(mB,Km,n,Input)
    
    return firstTerm + secondTerm

def repressor(mB,Km,n,Input):
    '''
        Repressor equation
    '''
    try:
        return mB/(1+(Input/Km)**n)
    #if Km is zero
    except ZeroDivisionError:
        try:
            return mB*(Km**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0
def activator(mB,Km,n,Input):
    '''
        Activator equation
    '''
    try:
        return mB/(1+(Km/Input)**n)
    #if Input is zero
    except ZeroDivisionError:
        try:
            return mB*(Input**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0
