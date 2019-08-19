#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

qdots = [.01,.05,.075,.1,.15,.2]

qdot = .05 # W
capthex = 1. # K
capti = capthex
rho = 145 # kg/m^3 -- density of He-II at ~1 K
v = 8./1000. # m^3 -- volume of He-II bottle
c = 100. # J/(kg*K) -- specific heat capacity
m = 3. # dimensionless -- the GM exponent

dx = 0.008 # m -- "length" of channel (hole)
ftinv = 377028198.7 # W3/(m5*K) -- conductivity function at 1 K
ahole = dx*dx # m^2 -- effective area of hole

k = ftinv*ahole**m/dx # W3/K

dt = 1. # s time step
nsteps = 1000

t = 0 # time
capt = capti # temperature

ts = np.empty(nsteps)
capts = np.empty(nsteps)

for i in np.arange(nsteps):
    ts[i]=t
    capts[i]=capt
    dcapt=dt*(qdot-((capt-capthex)*k)**(1/m))/(rho*v*c)
    t=t+dt
    capt=capt+dcapt
    print(i,t,dcapt,capt)

plt.plot(ts,capts)
plt.show()
