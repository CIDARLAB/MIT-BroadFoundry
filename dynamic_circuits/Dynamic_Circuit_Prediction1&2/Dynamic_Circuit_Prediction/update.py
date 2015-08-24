import inputs

#Change Functions
def getProteinChange(initProteinVal,initmRNAVal,gateProperties):
    #dX_dt = pB*mX - ap*Xi
    return gateProperties['pB']*initmRNAVal - gateProperties['ap']*initProteinVal
    
def getmRNAChange(t,initmRNAVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,initProtein):
    #degradation rate
    answer = -gateProperties['am']*initmRNAVal

    #production rate from each promoter
    fanIn = gateProperties['INPUT']
    for i in range(len(fanIn)):
        InputName = fanIn[i]
        indexOfInput = input_and_logic_gate_names.index(InputName)
        InputProperties = input_and_logic_gate_dictionaries[indexOfInput]
        mB = InputProperties['mB']
        Km = InputProperties['Km']
        n = InputProperties['n']

        #If this is another gate then just get the concentration of the protein from it.
        if InputProperties['TYPE'] != 'INPUT':
            Input = initProtein[logic_gate_names.index(InputName)]
        #If it is an input, then use the equation specified for the input to determine
        #the concentration of the input at that time step.
        elif InputProperties['TYPE'] == 'INPUT':
            Inputfunc = InputProperties['INPUT']
            Input = Inputfunc[0](t,*Inputfunc[1:])
        if InputProperties['INPUT_EFFECT']=="REPRESS":
            answer += repressor(mB,Km,n,Input)
        elif InputProperties['INPUT_EFFECT']=="ACTIVATE":
            answer += activator(mB,Km,n,Input)    
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