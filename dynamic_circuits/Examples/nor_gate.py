import matplotlib.pyplot as plt
import inputs
from scipy.integrate import odeint
import numpy as np

xMin = 0
#yMin = 0
xMax = 200
#yMax = 55

mXi = 0
mYi = 0
Xi = 0
Yi = 0
initc = [mXi, Xi, mYi, Yi]

M_HALFLIFE = 2.             # half-life of mRNA transcripts
P_HALFLIFE = 10             # half-life of protein

My_HALFLIFE = 2.             # half-life of mRNA transcripts
Py_HALFLIFE = 10             # half-life of protein

amx = np.log(2)/M_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
amy = np.log(2)/M_HALFLIFE  # degradation rate of mRNA (half-life = 2 min)
mBx = 30                    # induced rate of transcription (30 transcript/min)
mBy = 30                    # induced rate of transcription (30 transcript/min)

n = 2                       # hill coefficient
Km = 40                     # repression threshold of promoters, monomers per cell

itr = 500                  # time iterations
t = np.linspace(xMin, xMax, itr)   # time grid

ax = np.log(2)/P_HALFLIFE   # degradation rate of protein (half-life = 10 min)
ay = np.log(2)/P_HALFLIFE   # degradation rate of protein (half-life = 10 min)
Bx = 20.                    # translation efficiency (20 proteins/transcript)
By = 20.                    # translation efficiency (20 proteins/transcript)

# System of DEs and input equation.
def f(state, t):
        #concentration of mRNA and protein at this state
        mX = state[0]
        Xi = state[1]
        mY = state[2]
        Yi = state[3]
        
        #IPTG input        
        IPTGi = inputs.linInput(t, 0, 10)
        aTci = inputs.linInput(t, 0, 0)
        
        #mRNA production
        dmX_dt = -amx*mX + mBx*(np.power(IPTGi,n)/(np.power(IPTGi,n)+np.power(Km,n))+np.power(aTci,n)/(np.power(aTci,n)+np.power(Km,n)))
        dmY_dt = -amy*mY + mBy/(1+np.power(Xi/Km,n))
        
        #protein production
        dX_dt = Bx*mX - ax*Xi
        dY_dt = By*mY - ay*Yi
        
        return [dmX_dt, dX_dt, dmY_dt, dY_dt]
                
soln_f = odeint(f, initc, t)
#soln_g = odeint(g, initc, t)

Xm_f = soln_f[:,0]
Xp_f = soln_f[:,1]
Ym_f = soln_f[:,2]
Yp_f = soln_f[:,3]

plt.figure()
#plt.axis([xMin, xMax, yMin, yMax])
plt.plot(t,Xm_f, 'r--', label = 'Xm_f')
plt.plot(t,Ym_f, 'b--', label = 'Ym_f')
plt.xlabel('Time (min)')
plt.ylabel('Transcripts per cell')
plt.title('Example mRNA Data')
plt.legend(loc=0)

plt.figure()
#plt.axis([xMin, xMax, 0, 500])
plt.plot(t,Xp_f, 'r-', label = 'X_f')
plt.plot(t,Yp_f, 'b-', label = 'Y_f')
plt.xlabel('Time (min)')
plt.ylabel('Proteins per cell')
plt.title('Example Protein Data')
plt.legend(loc=0)