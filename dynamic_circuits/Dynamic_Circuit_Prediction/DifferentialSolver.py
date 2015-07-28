# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 18:19:43 2015

@author: Arinze
"""
import numpy as np

def differentialSolver(func,initc,timeArray,args=(),n2=None):
#    if n2==None:
#        n2=len(timeArray)/100
    results = np.asarray(initc)
    initc = np.asarray(initc)
    theRange = range(1,len(timeArray))
    for i in theRange:
        
        dt = timeArray[i]-timeArray[i-1]
        change = np.asarray(func(initc,timeArray[i],*args))
        for j in range(len(change)):
            initc = initc + change*dt
        results = np.vstack((results,initc))

#        if nLastTermsEqual(n2+1,results):
#            dt2 = timeArray[i+n2]-timeArray[i]
#            change2 = np.asarray(func(initc,timeArray[i+n2],*args))
#            change2 = np.around(change2,decimals=7)
#            initc2 = initc + change2*dt2
#            if (initc==initc2).all():
#                for j in range(1,n2+1):
#                    results = np.vstack((results,initc2))
#                    theRange.remove(j+i)
#                theRange.remove(j+1+i)
                    
        
        if i%1000 ==0:
            print 100*float(i)/len(timeArray),"%"
    return results
    
def nLastTermsEqual(n,array):
    arrayLength = len(array)
    if arrayLength<n:
        return False
    stock = array[-1]
    for i in range(arrayLength-n,arrayLength):
        if not (array[i]==stock).all():
            return False
    return True
        