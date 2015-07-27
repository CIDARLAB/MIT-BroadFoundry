import inputs

#Change Functions
def getProteinChange(initProteinVal,initmRNAVal,gateProperties):
    #dX_dt = pB*mX - ap*Xi
    return gateProperties['pB']*initmRNAVal - gateProperties['ap']*initProteinVal
    
def getmRNAChange(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    #first term is degradation rate other terms are from
    firstTerm = -gateProperties['am']*initmRNAVal
    
    #These gates only have one input so they only have two terms
    if gateProperties['TYPE']=='BUFFER' or gateProperties['TYPE']=='NOT':
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
        if gateProperties['INPUT_EFFECT']=="REPRESS":
            secondTerm = repressor(mB,Km,n,Input)
        elif gateProperties['INPUT_EFFECT']=="ACTIVATE":
            secondTerm = activator(mB,Km,n,Input)
        answer =  firstTerm + secondTerm   
        
    #These gates have two inputs so they will have a total of three terms
    elif gateProperties['TYPE']=='OR' or gateProperties['TYPE']=='NOR':
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
        if gateProperties['INPUT1_EFFECT']=="REPRESS":
            secondTerm = repressor(mB1,Km1,n1,Input1)
        elif gateProperties['INPUT1_EFFECT']=="ACTIVATE":
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
        if gateProperties['INPUT2_EFFECT']=="REPRESS":
            thirdTerm = repressor(mB2,Km2,n2,Input2)
        elif gateProperties['INPUT2_EFFECT']=="ACTIVATE":
            thirdTerm = activator(mB2,Km2,n2,Input2)
        answer = firstTerm + secondTerm + thirdTerm        
        
    return answer

#Repressor and Activator Hill Equations
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