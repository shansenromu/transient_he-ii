#!/usr/bin/python
from scipy.optimize import curve_fit
import numpy as np
#import matplotlib.pyplot as plt
#plt.style.use('seaborn-whitegrid')
import ROOT
from array import array
#import plot_heater_test_function as ht

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

pdf = "heat-bg-transeints.pdf"

canvas = ROOT.TCanvas('c', 'c')

#qdots = [0.,.010,.025,.050,.075,.100,.150,.200,.250] # W
qdots  =[0.,0.025,0.075,0.100,0.200]#W

qdot_bg = 0.198 # W

capthex = 0.89 # K

rho = 145. # kg/m^3 -- density of He-II at ~1 K
v = 8.5/1000. # m^3 -- volume of He-II bottle
m = 3. # dimensionless -- the GM exponent

dx = 0.0067 # m -- "length" of channel (hole)
dy = dx # m -- side length of hole area
ahole = dy*dy # m^2 -- effective area of hole


def FillGraphs(heat,step,graph1,graph2,graph3,graph4,graph5,time,cap):
    if heat == 0.:
        graph1.SetPoint(step,time,cap)
    if heat == 0.025:
        graph2.SetPoint(step,time,cap)
    if heat == 0.075:
        graph3.SetPoint(step,time,cap)
    if heat == 0.100:
        graph4.SetPoint(step,time,cap)
    if heat == 0.200:
        graph5.SetPoint(step,time,cap)

def func(t,a,b,c):
    return a*np.exp(-t/b)+c

def c(capt): # J/(kg*K)
    #return 100. # previous model
    if (capt<=.6):
        return 20.4*capt**3 # van Sciver Eq. (6.28)
    elif (capt<1.1):
        return 108.*capt**6.7 # van Sciver Eq. (6.29a)
    else:
        return 117.*capt**5.6 # van Sciver Eq. (6.29b)

def ftinv_sciver(capt):
    t_lambda=2.17 # K
    smallt=capt/t_lambda # dimensionless
    multiplier=(smallt**5.7*(1.-smallt**5.7))**3 # dimensionless
    s_lambda=1559. # J/(kg K)
    A_lambda=1450. # (m s)/kg
    rho=145. # kg/m^3
    g_lambda=rho**2*s_lambda**4*t_lambda**3/A_lambda # W^3/(K m^5)
    return g_lambda*multiplier


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

#fig1=plt.figure()
#Initialize the plots
gr_steps = nsteps
gr0mW = ROOT.TGraph(gr_steps,np.empty(gr_steps),np.empty(gr_steps))
gr25mW = ROOT.TGraph(gr_steps,np.empty(gr_steps),np.empty(gr_steps))
gr75mW = ROOT.TGraph(gr_steps,np.empty(gr_steps),np.empty(gr_steps))
gr100mW = ROOT.TGraph(gr_steps,np.empty(gr_steps),np.empty(gr_steps))
gr200mW = ROOT.TGraph(gr_steps,np.empty(gr_steps),np.empty(gr_steps))

for qdot in qdots:

    t = 0 # time
    capt = captstart # temperature

    ts = np.empty(nsteps)
    capts = np.empty(nsteps)

    ftinv_integral=dx*((qdot_bg)/ahole)**m
    for i in np.arange(nsteps):
        ts[i]=t
        capts[i]=capt
        FillGraphs(qdot,int(i),gr0mW,gr25mW,gr75mW,gr100mW,gr200mW,float(t),float(capt))
        dcapt=dt*(qdot+qdot_bg-ahole*(ftinv_integral/dx)**(1/m))/(rho*v*c(capt))
        t=t+dt
        capt=capt+dcapt
        ftinv_integral=ftinv_integral+dcapt*ftinv_sciver(capt)
        #print(i,t,dcapt,capt,ftinv_integral)

    #plt.plot(ts,capts,label='%4.3f W'%qdot)
    #intead of this I am filling above


    # # fit curve
    # popt,pcov=curve_fit(func,ts,capts)
    # print(qdot)
    # print(popt)
    # #This curve fit is unnessary at the moment
    # #plt.plot(ts,func(ts,*popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    #
    ###
    # Fit a piece-wise curve (since the taus would never be the same, and Delta T is under reported)
    ###
    



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

#plt.xlabel('Time (s)')
gr0mW.GetXaxis().SetTitle('Time (s)')
gr25mW.GetXaxis().SetTitle('Time (s)')
gr75mW.GetXaxis().SetTitle('Time (s)')
gr100mW.GetXaxis().SetTitle('Time (s)') 
gr200mW.GetXaxis().SetTitle('Time (s)') 
#plt.ylabel('Bottle temperature (K)')
gr0mW.GetYaxis().SetTitle('Bottle temperature (K)')
gr25mW.GetYaxis().SetTitle('Bottle temperature (K)')
gr75mW.GetYaxis().SetTitle('Bottle temperature (K)')
gr100mW.GetYaxis().SetTitle('Bottle temperature (K)')
gr200mW.GetYaxis().SetTitle('Bottle temperature (K)')
#plt.legend()
#LEGEND GOES HERE
gr0mW.Draw('al')
gr25mW.Draw('same')
gr75mW.Draw('same')
gr100mW.Draw('same')
gr200mW.Draw('same')
canvas.Print(pdf+'(')

#for qdot in qdots:
#fig2=plt.figure()
#plt.plot(qdots,captend,label='after 1000 s')
#plt.plot(qdots,captendthy,label='at infinity')
gr_captend = ROOT.TGraph(len(qdots),np.array(qdots),captend)
gr_captendthy = ROOT.TGraph(len(qdots),np.array(qdots),captendthy)
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('Power (W)')
gr_captend.GetXaxis().SetTitle('Power (W)')
gr_captendthy.GetXaxis().SetTitle('Power (W)')
#plt.ylabel('Final Temperature (K)')
gr_captend.GetXaxis().SetTitle('Final Temperature (K)')
gr_captendthy.GetXaxis().SetTitle('Final Temperature (K)')
#plt.legend()
#LEGEND GOES HERE
#plt.show()
gr_captend.Draw('al')
gr_captendthy.Draw('same')
canvas.Print(pdf+')')
