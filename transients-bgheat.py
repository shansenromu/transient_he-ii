#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

qdots = [.01,.05,.075,.1,.15,.2,.25] # W

qdot_bg = .15 # W

qdot = .05 # W
capthex = 1. # K
rho = 145. # kg/m^3 -- density of He-II at ~1 K
v = 8./1000. # m^3 -- volume of He-II bottle
c = 100. # J/(kg*K) -- specific heat capacity
m = 3. # dimensionless -- the GM exponent

def ftinv_sciver(capt):
    t_lambda=2.17 # K
    smallt=capt/t_lambda # dimensionless
    multiplier=(smallt**5.7*(1.-smallt**5.7))**3 # dimensionless
    s_lambda=1559. # J/(kg K)
    A_lambda=1450. # (m s)/kg
    rho=145. # kg/m^3
    g_lambda=rho**2*s_lambda**4*t_lambda**3/A_lambda # W^3/(K m^5)
    return g_lambda*multiplier

dx = 0.01 # m -- "length" of channel (hole)
ahole = dx*dx # m^2 -- effective area of hole


# First, calculate starting temperature by assuming it has equilibrated.

ftinv_integral=0
ftinv_goal=dx*(qdot_bg/ahole)**3
captprime=capthex
deltacaptstep=0.0001
while (ftinv_integral<ftinv_goal):
    ftinv_integral=ftinv_integral+deltacaptstep*ftinv_sciver(captprime)
    captprime=captprime+deltacaptstep
captstart=captprime
print('Starting temperature is %f'%captstart)

dt = .01 # s time step
nsteps = 100000

iq=0
captend=np.empty(len(qdots))
captendthy=np.empty(len(qdots))
fig1=plt.figure()
for qdot in qdots:

    t = 0 # time
    capt = captstart # temperature

    ts = np.empty(nsteps)
    capts = np.empty(nsteps)

    ftinv_integral=0.
    for i in np.arange(nsteps):
        ts[i]=t
        capts[i]=capt
        dcapt=dt*(qdot+qdot_bg-ahole*(ftinv_integral/dx)**(1/m))/(rho*v*c)
        t=t+dt
        capt=capt+dcapt
        ftinv_integral=ftinv_integral+dcapt*ftinv_sciver(capt)
        #print(i,t,dcapt,capt,ftinv_integral)

    plt.plot(ts,capts,label='%4.3f W'%qdot)

    captend[iq]=capt

    # Calculate theoretical balance of qdot and ftinv integral

    ftinv_integral=0
    ftinv_goal=dx*((qdot+qdot_bg)/ahole)**3
    captprime=captstart
    deltacaptstep=0.001
    while (ftinv_integral<ftinv_goal):
        ftinv_integral=ftinv_integral+deltacaptstep*ftinv_sciver(captprime)
        captprime=captprime+deltacaptstep

    captendthy[iq]=captprime
    iq=iq+1

plt.xlabel('Time (s)')
plt.ylabel('Bottle temperature (K)')
plt.legend()

#for qdot in qdots:
fig2=plt.figure()
plt.plot(qdots,captend,label='after 1000 s')
plt.plot(qdots,captendthy,label='at infinity')
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Power (W)')
plt.ylabel('Final Temperature (K)')
plt.legend()
plt.show()
