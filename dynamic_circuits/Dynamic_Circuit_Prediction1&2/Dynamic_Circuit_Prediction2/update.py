
#Change Functions
def getREUChange(t,initREUVal,gateProperties,input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,REUVals):
    #first term is degradation rate other terms are from
    degRate = gateProperties['degRate']
    
    
    #equation: dY/dt = Bo + B/(1+(x/K)^n) - alpha*Y
    #where: Bo = alpha*Pmin
    #       B = alpha*(Pmax-Pmin)
    
    
    ### minus alpha*Y    
    answer = -degRate*initREUVal

    #REUVals has the current value for each protein    
    
    
    #add together the effects of each input

    #fanIn is an array of gate names.
    fanIn = gateProperties['INPUT']
    
    for i in range(len(fanIn)):
        
        #input to the gate (not circuit input, necessarily)
        InputName = fanIn[i]
        
        indexOfInput = input_and_logic_gate_names.index(InputName)
        InputProperties = input_and_logic_gate_dictionaries[indexOfInput]
        Pmin = InputProperties['Pmin']
        Pmax = InputProperties['Pmax']
        Km = InputProperties['Km']
        n = InputProperties['n']

        #Bo.
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
            #hill fn
            answer += degRate*repressor(Pmax-Pmin,Km,n,Input)
        elif InputProperties['INPUT_EFFECT']=="ACTIVATE":
            #hill fn
            answer += degRate*activator(Pmax-Pmin,Km,n,Input) 
        
        #return array of [dX/dt, dY/dt, etc.]
    return answer
    
#Repressor and Activator Hill Equations
def repressor(mB,Km,n,Input):
    '''
        Repressor equation
    '''
    #Prevent imaginary numbers. Also, Input should never be below zero.
    if Input<0:
        Input = 0
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
    #Prevent imaginary numbers. Also, Input should never be below zero.
    if Input<0:
        Input = 0
    try:
        return mB/(1+(Km/Input)**n)
    #if Input is zero
    except ZeroDivisionError:
        try:
            return mB*(Input**n)/(Input**n+Km**n)
        #if both Km and Input are zero
        except ZeroDivisionError:
            return 0