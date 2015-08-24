import inputs

#Change Functions
def getREUChange(t,initREUVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,REUVals):
    #first term is degradation rate other terms are from
    degRate = gateProperties['degRate']
    answer = -degRate*initREUVal
    
    #add together the effects of each input
    fanIn = gateProperties['INPUT']
    for i in range(len(fanIn)):
        InputName = fanIn[i]
        indexOfInput = input_and_logic_gate_names.index(InputName)
        InputProperties = input_and_logic_gate_dictionaries[indexOfInput]
        Pmin = InputProperties['Pmin']
        Pmax = InputProperties['Pmax']
        Km = InputProperties['Km']
        n = InputProperties['n']

        answer += Pmin*degRate
        #If this is another gate then just get the concentration of the protein from it.
        if InputProperties['TYPE'] != 'INPUT':
            Input = REUVals[logic_gate_names.index(InputName)]
        #If it is an input, then use the equation specified for the input to determine
        #the concentration of the input at that time step.
        elif InputProperties['TYPE'] == 'INPUT':
            Inputfunc = InputProperties['INPUT']
            Input = Inputfunc[0](t,*Inputfunc[1:])
        if InputProperties['INPUT_EFFECT']=="REPRESS":
            answer += degRate*repressor(Pmax-Pmin,Km,n,Input)
        elif InputProperties['INPUT_EFFECT']=="ACTIVATE":
            answer += degRate*activator(Pmax-Pmin,Km,n,Input) 
        
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