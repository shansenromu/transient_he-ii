import ROOT
import math



def heater_test_function( x, p ):
  # x[0] x value (time axis)
  # p[0] Delta t1 amplitude of temp jump
  # p[1] tau1     rising time constant (s)
  # p[2] m        background slope (K/s) added relative to turn on time
  # p[3] b        background temperature (K)
  # p[4] t1       turn on time (s)
  # p[5] t2       turn off time (s)
  # p[6] Delta t2 amplitude of temp fall
  # p[7] tau2     falling time constant (s)
  
  if x[0] < p[4]:  # before heater turns on
    return p[3] + p[2]*(x[0]-p[4]);
  elif x[0] < p[5]: #while heater is on
    retval = p[3] + p[2]*(x[0]-p[4]);
    retval = retval + p[0]*(1 - math.exp( -(x[0]-p[4])/p[1] ) )
    return retval
  else:
    # heater turns off
    bprime = p[3] + p[0]*( 1 - math.exp( -(p[5]-p[4])/p[1] ) ) - p[6]
    retval = bprime + p[2]*(x[0]-p[4]);
    retval = retval + p[6]*math.exp( -(x[0]-p[5])/p[7] )
    return retval


def htf_sameEXP( x, p ): #heater test function with same exponetial values
  # x[0] x value (time axis)
  # p[0] Delta t1 amplitude of temp jump
  # p[1] tau1     rising time constant (s)
  # p[2] m        background slope (K/s) added relative to turn on time
  # p[3] b        background temperature (K)
  # p[4] t1       turn on time (s)
  # p[5] t2       turn off time (s)
  # p[6] Delta t2 amlitude of temp fall
  
  if x[0] < p[4]:  # before heater turns on
    return p[3] + p[2]*(x[0]-p[4]);
  elif x[0] < p[5]: #while heater is on
    retval = p[3] + p[2]*(x[0]-p[4]);
    retval = retval + p[0]*(1 - math.exp( -(x[0]-p[4])/p[1] ) )
    return retval
  else:
    # heater turns off
    bprime = p[3] + p[0]*( 1 - math.exp( -(p[5]-p[4])/p[1] ) ) - p[6]
    retval = bprime + p[2]*(x[0]-p[4]);
    retval = retval + p[6]*math.exp( -(x[0]-p[5])/p[1] )
    return retval

def htf_sameDT( x, p ): #heater test function with same exponetial values
  # x[0] x value (time axis)
  # p[0] Delta t1 amplitude of temp jump
  # p[1] tau1     rising time constant (s)
  # p[2] m        background slope (K/s) added relative to turn on time
  # p[3] b        background temperature (K)
  # p[4] t1       turn on time (s)
  # p[5] t2       turn off time (s)
  #               Delta t2 amlitude of temp fall
  # p[6] tau2
  if x[0] < p[4]:  # before heater turns on
    return p[3] + p[2]*(x[0]-p[4]);
  elif x[0] < p[5]: #while heater is on
    retval = p[3] + p[2]*(x[0]-p[4]);
    retval = retval + p[0]*(1 - math.exp( -(x[0]-p[4])/p[1] ) )
    return retval
  else:
    # heater turns off
    bprime = p[3] + p[0]*( 1 - math.exp( -(p[5]-p[4])/p[1] ) ) - p[0]
    retval = bprime + p[2]*(x[0]-p[4]);
    retval = retval + p[0]*math.exp( -(x[0]-p[5])/p[6] )
    return retval

def htf_sameTexp( x, p ): #heater test function with same exponetial values
  # x[0] x value (time axis)
  # p[0] Delta t1 amplitude of temp jump
  # p[1] tau1     rising time constant (s)
  # p[2] m        background slope (K/s) added relative to turn on time
  # p[3] b        background temperature (K)
  # p[4] t1       turn on time (s)
  # p[5] t2       turn off time (s)
  #               Delta t2 amlitude of temp fall
  #               tau2=tau1
  
  if x[0] < p[4]:  # before heater turns on
    return p[3] + p[2]*(x[0]-p[4]);
  elif x[0] < p[5]: #while heater is on
    retval = p[3] + p[2]*(x[0]-p[4]);
    retval = retval + p[0]*(1 - math.exp( -(x[0]-p[4])/p[1] ) )
    return retval
  else:
    # heater turns off
    bprime = p[3] + p[0]*( 1 - math.exp( -(p[5]-p[4])/p[1] ) ) - p[0]
    retval = bprime + p[2]*(x[0]-p[4]);
    retval = retval + p[0]*math.exp( -(x[0]-p[5])/p[1] )
    return retval



def cool_function(x, p): #cooling portion only
    # x[0] x value (time axis)
    # p[0] Delta t1 amplitude of temp fall
    # p[1] tau1     rising time constant (s)
    # p[2] m        background slope (K/s) added relative to turn on time
    # p[3] b        background temperature (K)
    # p[4] t1       turn off time (s)


    if x[0] < p[4]: # while heater turns on
        return p[3] + p[2]*(x[0]-p[4])+p[0];
    else:
        # heater turns off
        retval = p[3] + p[2]*(x[0]-p[4]);
        retval = retval + p[0]*math.exp( -(x[0]-p[4])/p[1] )
        return retval

def high_function( x, p ):
  # x[0] x value (time axis)
  # p[0] Delta t1 amplitude of temp jump
  # p[1] tau1     rising time constant (s)
  #### p[2] m        background slope (K/s) added relative to turn on time
  # p[2] b        background temperature (K)
  # p[3] t1       turn on time (s)

  
  if x[0] < p[3]:  # before heater turns on
    return p[2]# + p[2]*(x[0]-p[4])
  else: #while heater is on
    retval = p[2] #+ p[2]*(x[0]-p[4]);
    retval = retval + p[0]*(1 - math.exp( -(x[0]-p[3])/p[1] ) )
    return retval



def linear_function( x, p ): #heater test function with same exponetial values
  # x[0] x value (time axis)
  # p[0] m        background slope (K/s) added relative to turn on time
  # p[1] b        background temperature (K)
  return p[1] + p[0]*(x[0]);



def SingleExpoWithSlopeBackground(oS, pp0 ,pp1 , pp2, pp3):
    SingleExpoWithSlopeBackground = ROOT.TF1('SingleExpoWithBackground', '[0]*exp((-x+{0})/[1]) + [2]*(x-{0})+[3]'.format(oS))
    SingleExpoWithSlopeBackground.SetParameters(pp0,pp1,pp2,pp3 )
    SingleExpoWithSlopeBackground.SetParName(0, '#Delta Temp')
    SingleExpoWithSlopeBackground.SetParName(1, 'time constant #tau')
    SingleExpoWithSlopeBackground.SetParName(2, 'Background Slope')
    SingleExpoWithSlopeBackground.SetParName(3, 'Final Temp')
    SingleExpoWithSlopeBackground.SetParLimits(0, -1, 1)
    SingleExpoWithSlopeBackground.SetParLimits(1, -1000, 1000)
    SingleExpoWithSlopeBackground.SetParLimits(2, -1, 1)
    SingleExpoWithSlopeBackground.SetParLimits(3, 0, 2)
    return SingleExpoWithSlopeBackground



def HeaterTestFunction(pars, tmin, tmax):
    HeaterTestFunction = ROOT.TF1('HeaterTest', heater_test_function, tmin, tmax, 8 )
    HeaterTestFunction.SetParameters( pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7] )
    HeaterTestFunction.SetParName(0, '#Delta T_{1} (K)')
    HeaterTestFunction.SetParName(1, '#tau_{1} (s)    ')
    HeaterTestFunction.SetParName(2, 'BG slope (K/s)  ')
    HeaterTestFunction.SetParName(3, 'Offset T        ')
    HeaterTestFunction.SetParName(4, 'turn on (s)     ')
    HeaterTestFunction.SetParName(5, 'turn off (s)    ')
    HeaterTestFunction.SetParName(6, '#Delta T_{2} (K)')
    HeaterTestFunction.SetParName(7, '#tau_{2} (s)    ')
    HeaterTestFunction.SetParLimits(0, 0, 1.0)
    HeaterTestFunction.SetParLimits(1, 1, 1000)
    HeaterTestFunction.SetParLimits(2, -0.1, 0.1)
    HeaterTestFunction.SetParLimits(3, 0, 2)    
    HeaterTestFunction.SetParLimits(4, 0, 14000)    
    HeaterTestFunction.SetParLimits(5, 0, 14000)    
    HeaterTestFunction.SetParLimits(6, 0, 1.0)
    HeaterTestFunction.SetParLimits(7, 1, 1000)
    HeaterTestFunction.SetNpx(300)
    return HeaterTestFunction

def HTF_Same_exp(pars, tmin, tmax):
    HTF_Same_exp = ROOT.TF1('HeaterTest', htf_sameEXP, tmin, tmax, 7 )
    HTF_Same_exp.SetParameters( pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6] )

    HTF_Same_exp.SetParName(0, '#Delta T_{1} (K)')
    HTF_Same_exp.SetParName(1, '#tau_{1} (s)    ')
    HTF_Same_exp.SetParName(2, 'BG slope (K/s)  ')
    HTF_Same_exp.SetParName(3, 'Offset T        ')
    HTF_Same_exp.SetParName(4, 'turn on (s)     ')
    HTF_Same_exp.SetParName(5, 'turn off (s)    ')
    HTF_Same_exp.SetParName(6, '#Delta T_{2} (K)')
    HTF_Same_exp.SetParLimits(0, 0, 1.0)
    HTF_Same_exp.SetParLimits(1, 1, 1000)
    HTF_Same_exp.SetParLimits(2, -0.1, 0.1)
    HTF_Same_exp.SetParLimits(3, 0, 2)    
    HTF_Same_exp.SetParLimits(4, 0, 14000)    
    HTF_Same_exp.SetParLimits(5, 0, 14000)    
    HTF_Same_exp.SetParLimits(6, 0, 1.0)
    HTF_Same_exp.SetNpx(300)
    return HTF_Same_exp

def HTF_Same_T(pars, tmin, tmax):
    HTF_Same_exp = ROOT.TF1('HeaterTest', htf_sameDT, tmin, tmax, 7 )
    HTF_Same_exp.SetParameters( pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6] )

    HTF_Same_exp.SetParName(0, '#Delta T_{1} (K)')
    HTF_Same_exp.SetParName(1, '#tau_{1} (s)    ')
    HTF_Same_exp.SetParName(2, 'BG slope (K/s)  ')
    HTF_Same_exp.SetParName(3, 'Offset T        ')
    HTF_Same_exp.SetParName(4, 'turn on (s)     ')
    HTF_Same_exp.SetParName(5, 'turn off (s)    ')
    HTF_Same_exp.SetParName(6, '#tau_{2} (s)    ')
    HTF_Same_exp.SetParLimits(0, 0, 1.0)
    HTF_Same_exp.SetParLimits(1, 1, 1000)
    HTF_Same_exp.SetParLimits(2, -0.1, 0.1)
    HTF_Same_exp.SetParLimits(3, 0, 2)    
    HTF_Same_exp.SetParLimits(4, 0, 14000)    
    HTF_Same_exp.SetParLimits(5, 0, 14000)    
    HTF_Same_exp.SetParLimits(6, 0, 1000)
    HTF_Same_exp.SetNpx(300)
    return HTF_Same_exp

def HTF_Same_Texp(pars, tmin, tmax):
    HTF_Same_exp = ROOT.TF1('HeaterTest', htf_sameTexp, tmin, tmax, 6 )
    HTF_Same_exp.SetParameters( pars[0], pars[1], pars[2], pars[3], pars[4], pars[5])

    HTF_Same_exp.SetParName(0, '#Delta T_{1} (K)')
    HTF_Same_exp.SetParName(1, '#tau_{1} (s)    ')
    HTF_Same_exp.SetParName(2, 'BG slope (K/s)  ')
    HTF_Same_exp.SetParName(3, 'Offset T        ')
    HTF_Same_exp.SetParName(4, 'turn on (s)     ')
    HTF_Same_exp.SetParName(5, 'turn off (s)    ')
    HTF_Same_exp.SetParLimits(0, 0, 1.0)
    HTF_Same_exp.SetParLimits(1, 1, 1000)
    HTF_Same_exp.SetParLimits(2, -0.1, 0.1)
    HTF_Same_exp.SetParLimits(3, 0, 2)    
    HTF_Same_exp.SetParLimits(4, 0, 14000)    
    HTF_Same_exp.SetParLimits(5, 0, 14000)    
    HTF_Same_exp.SetNpx(300)
    return HTF_Same_exp


def Cool_Test(pars, tmin, tmax):
    Cool_Test = ROOT.TF1('HeaterTest', cool_function, tmin, tmax, 5 )
    Cool_Test.SetParameters( pars[0], pars[1], pars[2], pars[3], pars[4] )
    Cool_Test.SetParName(0, '#Delta T_{2} (K)')
    Cool_Test.SetParName(1, '#tau_{2} (s)    ')
    Cool_Test.SetParName(2, 'BG slope (K/s)  ')
    Cool_Test.SetParName(3, 'Offset T        ')
    Cool_Test.SetParName(4, 'turn off (s)    ')
    Cool_Test.SetParLimits(0, 0, 1.0)
    Cool_Test.SetParLimits(1, 1, 1000)
    Cool_Test.SetParLimits(2, -0.1, 0.1) 
    Cool_Test.SetParLimits(3, 0, 2)    
    Cool_Test.SetParLimits(4, 0, 14000)
    Cool_Test.SetNpx(300)
    return Cool_Test

def High_Test(pars, tmin, tmax):
    High_Test = ROOT.TF1('HeaterTest', high_function, tmin, tmax, 4 )
    High_Test.SetParameters( pars[0], pars[1], pars[2], pars[3])
    High_Test.SetParName(0, '#Delta T_{1} (K) ')
    High_Test.SetParName(1, '#tau_{1} (s)     ')
    #High_Test.SetParName(2, 'BG slope (K/s)   ')
    High_Test.SetParName(2, 'Starting Temp (K)')
    High_Test.SetParName(3, 'turn on (s)      ')
    High_Test.SetParLimits(0, 0, 1.0)
    High_Test.SetParLimits(1, 1, 1000)
    #High_Test.SetParLimits(2, -0.1, 0.1) 
    High_Test.SetParLimits(2, 0., 2.)    
    High_Test.SetParLimits(3, -10, 10)
    High_Test.SetNpx(300)
    return High_Test


def Linear_Function(pars, tmin, tmax):
    Linear_Test = ROOT.TF1('HeaterTest', linear_function, tmin, tmax, 2 )
    Linear_Test.SetParameters( pars[0], pars[1] )
    Linear_Test.SetParName(0, 'BG slope (K/s)  ')
    Linear_Test.SetParName(1, 'Offset T        ')
    Linear_Test.SetParLimits(0, -0.1, 0.1)  
    Linear_Test.SetParLimits(1, 0, 2)    
    return Linear_Test

#test1 = HeaterTestFunction( [0.024, 32.0, 5.0e-6, 1.16, 900.0, 1000.0, 0.024, 32.0 ], 600.0, 1200.0 )

#tc = ROOT.TCanvas()
#test1.Draw()
#tc.Print("plot_heater_test_function.png")






