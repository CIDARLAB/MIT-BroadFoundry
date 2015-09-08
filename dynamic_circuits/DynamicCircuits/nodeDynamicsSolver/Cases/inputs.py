"""
Created on Thu Aug 20 1:28:30 2015

@author: Alex Lim
"""

# t: integer value of the time point at which
# signal: a list of floats signifying the signal over time
# maxTime: integer value of the max time
def getSquareWavePoint(t, signal, maxTime):
    """
    Returns the value corresponding with the square wave at a given time
    """
    divide = float(maxTime)/len(signal)
    		
    if t < 0:
        return signal[0]
    elif int(t/divide) >= len(signal):
        return signal[-1]
    else:
        j = int(t/divide)
        return signal[j]