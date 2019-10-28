#!/usr/bin/python
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

qdots = [0.,.025,.05,.075,.1,.15,.2,.25] # W

qdot_bg = .15 # W

capthex = .9 # K
rho = 145. # kg/m^3 -- density of He-II at ~1 K
v = 8./1000. # m^3 -- volume of He-II bottle
m = 3. # dimensionless -- the GM exponent

def func(t,a,b,c):
    return a*np.exp(-t/b)+c

def c(t): # J/(kg*K)
    #return 100. # previous model
    if (t<=.6):
        return 20.4*t**3 # van Sciver Eq. (6.28)
    elif (t<1.1):
        return 108.*t**6.7 # van Sciver Eq. (6.29a)
    else:
        return 117.*t**5.6 # van Sciver Eq. (6.29b)

def ftinv_sciver(capt):
    t_lambda=2.17 # K
    smallt=capt/t_lambda # dimensionless
    multiplier=(smallt**5.7*(1.-smallt**5.7))**3 # dimensionless
    s_lambda=1559. # J/(kg K)
    A_lambda=1450. # (m s)/kg
    rho=145. # kg/m^3
    g_lambda=rho**2*s_lambda**4*t_lambda**3/A_lambda # W^3/(K m^5)
    return g_lambda*multiplier

dx = 0.005 # m -- "length" of channel (hole)
dy = 0.005 # m -- side length of hole area
ahole = dy*dy # m^2 -- effective area of hole


# First, calculate starting temperature by assuming it has equilibrated.

ftinv_integral=0.
ftinv_goal=dx*(qdot_bg/ahole)**m
captprime=capthex
deltacaptstep=0.0001
while (ftinv_integral<ftinv_goal):
    ftinv_integral=ftinv_integral+deltacaptstep*ftinv_sciver(captprime)
    captprime=captprime+deltacaptstep
captstart=captprime
print('Starting temperature is %f'%captstart)
print(ftinv_goal)

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

    ftinv_integral=dx*((qdot_bg)/ahole)**m
    for i in np.arange(nsteps):
        ts[i]=t
        capts[i]=capt
        dcapt=dt*(qdot+qdot_bg-ahole*(ftinv_integral/dx)**(1/m))/(rho*v*c(capt))
        t=t+dt
        capt=capt+dcapt
        ftinv_integral=ftinv_integral+dcapt*ftinv_sciver(capt)
        #print(i,t,dcapt,capt,ftinv_integral)

    plt.plot(ts,capts,label='%4.3f W'%qdot)

    # fit curve
    popt,pcov=curve_fit(func,ts,capts)
    print(qdot)
    print(popt)
    plt.plot(ts,func(ts,*popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    
    captend[iq]=capt

    # Calculate theoretical balance of qdot and ftinv integral

    ftinv_integral=dx*((qdot_bg)/ahole)**m
    ftinv_goal=dx*((qdot+qdot_bg)/ahole)**m
    captprime=captstart
    deltacaptstep=0.0001
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
