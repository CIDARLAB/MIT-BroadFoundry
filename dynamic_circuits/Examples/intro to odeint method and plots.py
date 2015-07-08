'''
Example of use of scipy.integrate odeint method
'''
# zombie apocalypse modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.ion()

'''
Initialize constants
'''
P = 0       # birth rate
d = 0.0001  # natural death percent (per day)
B = 0.0095  # transmission percent  (per day)
G = 0.0001  # resurrect percent (per day)
A = 0.0001  # destroy percent  (per day)

'''
Initialize timescale and initial condition
'''
# initial conditions
t  = np.linspace(0, 6., 1000)   # time grid
S0 = 500.               # initial population of living
Z0 = 0                  # initial population of zombies
R0 = 0                  # initial population of dead
y0 = [S0, Z0, R0]       # initial condition vector

'''
Create differential equations representing relations between variable
'''
# solve the system dy/dt = f(y, t)
def f(y, t):
        Si = y[0]
        Zi = y[1]
        Ri = y[2]
        # the model equations (see Munz et al. 2009)
        #d(living)/dt = Birth - zombie transmit virus to living - natural death of living
        f0 = P - B*Si*Zi - d*Si
        #d(zombie)/dt = zombie transmit virus to living + resurrected dead as zombie - living rekill zombie 
        f1 = B*Si*Zi + G*Ri - A*Si*Zi
        #d(dead)/dt = natural death of living + living rekill zombie - resurrected dead as zombie
        f2 = d*Si + A*Si*Zi - G*Ri
        return [f0, f1, f2]

'''
The odeint method solves groups of differential equations and returns the 
solutions of X, Y, and Z as a list of 1000 points to be plotted between 
0 and 6.
'''
# solve the DEs
soln = odeint(f, y0, t)
S = soln[:, 0]
Z = soln[:, 1]
R = soln[:, 2]

plt.figure()
plt.plot(t, S, label='Living')
plt.plot(t, Z, label='Zombies')
plt.xlabel('Days from outbreak')
plt.ylabel('Population')
plt.title('Zombie Apocalypse - No Daily Births')
plt.legend(loc=0)

# change initial conditions
R0 = 0.01*S0   # 1% of initial pop is dead
S0 = S0-R0
P  = 10        # 10 new births daily
y0 = [S0, Z0, R0]

# solve the DEs
soln = odeint(f, y0, t)
S = soln[:, 0]
Z = soln[:, 1]
R = soln[:, 2]

plt.figure()
plt.plot(t, S, label='Living')
plt.plot(t, Z, label='Zombies')
plt.xlabel('Days from outbreak')
plt.ylabel('Population')
plt.title('Zombie Apocalypse - 1% Init. Pop. is Dead; 10 Daily Births')
plt.legend(loc=0)
