#import brian_no_units

from brian import *
import scipy.io as io
import cPickle
import sys
from StringIO import StringIO
from scipy import sparse
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import pickle



#ring model connectivity
def conn(k,sig,C):
  if abs(k)<=0.5:
    value=C*exp(-0.5*(k/sig)**2)/sqrt(2*pi)/sig
  else:
    value=C*exp(-0.5*((abs(k)-1)/sig)**2)/sqrt(2*pi)/sig
  return value


### 
def simulation(loadconnections=True, name_conections='connections_sp_20000.npz', 
  saveconnections=False, save_name='connections_sp_20000',
  timesimulation=10, dt_clock=0.1, stim_on=2000*ms, stim_off=3000*ms, pos_stim=0.5,
  epsE=60, epsI=0, 
  N=20000, prop_e=0.8, prop_i=0.2,  
  tauE = 20*ms , tauI=10*ms, taua=3*ms, taun=50*ms, taug=4*ms, taud =200*ms, tauf=450*ms, 
  Ustp = 0.04, Ustp0 = 0.03, R0=20.0/second, Vt=20*mV, Vr=-3.33*mV, refE=0.*ms, refI=0.*ms, 
  gEEA=533.3*mV*ms, gEEN=533.3*mV*ms, gEIA=67.2*mV*ms, gEIN=7.4*mV*ms, gIE=-138.6*mV*ms, gII=-90.6*mV*ms, 
  sigmaEE=15, sigmaEI=17.5, sigmaIE=15, sigmaII=15, 
  extE=0.*mV, extI=0.*mV, 
  stimE=0.24*mV,  stimI = 0.*mV):
  ##############  
  K=int(N/40)
  ############## Simulation of spiking neurons
  ##############    
  ### Taus for the spiking rate model of persistent activity
  ### taua AMPA synapse decay
  ### taun NMDA synapse decay
  ### taug GABA synapse decay
  ### taud  # synaptic depression tau
  ### tauf # synaptic facilitation tau
  ### Vt spike threshold
  ### Vr  reset value
  ### refE  refractory periods
  ### refI refractory periods
  ### Conductances
  ### gEEA 
  ### gEEN
  ### gEIA
  ### gEIN
  ### gIE
  ### gII
  #############
  defaultclock.reinit()
  defaultclock.dt = dt_clock*ms
  runtime=timesimulation*second
  ### Equations STP
  ## STP equations (http://www.scholarpedia.org/article/Short-term_synaptic_plasticity)
  u_ss = Ustp*(1.+R0*tauf)/(1.+Ustp*R0*tauf)
  u_ss0= Ustp0*(1.+R0*tauf)/(1.+Ustp0*R0*tauf)
  synfact_ss = u_ss/(1.+taud*R0*u_ss)
  synfact_ss0 = u_ss0/(1.+taud*R0*u_ss0)
  ntaud=(u_ss-synfact_ss0)/(synfact_ss0*R0*u_ss) #to avoid scaling
  #these are intermediate calculations needed for the equations below
  NE=int(ceil(prop_e*N)) # number of excitatory neurons
  NI=int(floor(prop_i*N)) # number of inhibitory neurons
  KE=int(ceil(prop_e*K)) # number of excitatory inputs
  KI=int(floor(prop_i*K)) # number of inhibitory inputs
  sigEE=sigmaEE/360.0
  sigEI=sigmaEI/360.0
  sigIE=sigmaIE/360.0
  sigII=sigmaII/360.0
  epsE = epsE/360.0
  gEEA=gEEA/sqrt(KE)/taua*synfact_ss0/synfact_ss #*np.sqrt(0.03/Ustp)
  gEEN=gEEN/sqrt(KE)/taun*synfact_ss0/synfact_ss #*np.sqrt(0.03/Ustp)
  gEIA=gEIA/sqrt(KE)/taua
  gEIN=gEIN/sqrt(KE)/taun
  gIE=gIE/sqrt(KI)/taug
  gII=gII/sqrt(KI)/taug
  stimE=stimE*sqrt(KE)
  stimI=stimI*sqrt(KE)
  #equations for each neuron. xpre and satura are shadow variables that track the saturated state
  # of a neuron's outgoing NMDA synapses, and this is used as a "synaptic depression" effect on x
  # (see modulation='satura' below) so as to mimic the NMDA saturation that is modeled in 
  # Compte et al., 2000
  eqs = '''
  dV/dt = (gea+gen+gi-V+Ix+Iext)/tau : mV
  dgea/dt = -gea/taua          : mV
  dgen/dt = -gen/taun   : mV
  dgi/dt = -gi/taug             : mV
  tau : second
  Ix : mV
  Iext : mV
  '''
  ##
  networkE=NeuronGroup(NE,model=eqs,threshold=Vt,reset=Vr, refractory=refE)
  networkI=NeuronGroup(NI,model=eqs,threshold=Vt,reset=Vr, refractory=refI)
  networkE.tau=tauE
  networkI.tau=tauI
  networkE.Ix=1.66*sqrt(KE)*mV ##maintain the same S/N ratio (KE number conection per neuron (SPARSE))
  networkI.Ix=1.85*0.83*sqrt(KE)*mV
  networkE.V = Vr + rand(NE)*(Vt - Vr) #Vt-2.0*mV + rand(NE) * 2.0*mV
  networkI.V = Vr + rand(NI)*(Vt - Vr) #Vt-2.0*mV + rand(NI) * 2.0*mV
  networkE.Iext=extE #### External input (IE0 en el rate model)
  networkI.Iext=extI
  ##
  if loadconnections: ## quicker
    loader = np.load('/home/david/Desktop/brian_simulations/' + str(name_conections), allow_pickle=True)
    WW = csr_matrix((loader['CEEd'], loader['CEEi'], loader['CEEp']), shape=loader['CEEs'])
    WW[WW != 0] = 1
    C1=Connection(networkE, networkE, 'gea', weight=gEEA*WW)
    C2=Connection(networkE, networkE, 'gen', weight=gEEN*WW)
    WW = csr_matrix((loader['CEId'], loader['CEIi'], loader['CEIp']), shape=loader['CEIs'])
    WW[WW != 0] = 1
    C3=Connection(networkE, networkI, 'gea', weight=gEIA*WW)
    C4=Connection(networkE, networkI, 'gen', weight=gEIN*WW)
    WW = csr_matrix((loader['CIEd'], loader['CIEi'], loader['CIEp']), shape=loader['CIEs'])
    WW[WW != 0] = 1
    C5=Connection(networkI, networkE, 'gi', weight=gIE*WW)
    WW = csr_matrix((loader['CIId'], loader['CIIi'], loader['CIIp']), shape=loader['CIIs'])
    WW[WW != 0] = 1
    C6=Connection(networkI, networkI, 'gi', weight=gII*WW)
  else:
    fE = float(KE)/float(NE)
    fI = float(KI)/float(NI)
    C1=Connection(networkE, networkE, 'gea', weight=gEEA, sparseness=lambda i,j: conn(float(i)/NE-float(j)/NE,sigEE,fE))
    C2=Connection(networkE, networkE, 'gen', weight=gEEN/gEEA*C1.W)
    C3=Connection(networkE, networkI, 'gea', weight=gEIA, sparseness=lambda i,j: conn(float(i)/NE-float(j)/NI,sigEI,fE))
    C4=Connection(networkE, networkI, 'gen', weight=gEIN/gEIA*C3.W)
    C5=Connection(networkI, networkE, 'gi', weight=gIE, sparseness=lambda i,j: conn(float(i)/NI-float(j)/NE,sigIE,fI))
    C6=Connection(networkI, networkI, 'gi', weight=gII, sparseness=lambda i,j: conn(float(i)/NI-float(j)/NI,sigII,fI))
  ###########
  ###########
  ########### 
  if saveconnections: ## do this the if you run a simulation with different number of neurons for the first time
    CEE = csr_matrix(C1.W)
    CEI = csr_matrix(C3.W)
    CIE = csr_matrix(C5.W)
    CII = csr_matrix(C6.W)
    np.savez_compressed('/home/david/Desktop/brian_simulations/'+str(save_name), CEEd=CEE.data, CEEi=CEE.indices, CEEp=CEE.indptr, CEEs=CEE.shape,
               CEId=CEI.data, CEIi=CEI.indices, CEIp=CEI.indptr, CEIs=CEI.shape,
               CIEd=CIE.data, CIEi=CIE.indices, CIEp=CIE.indptr, CIEs=CIE.shape,
               CIId=CII.data, CIIi=CII.indices, CIIp=CII.indptr, CIIs=CII.shape)
  #########
  #########
  mystpEEA=STP(C1,taud=taud,tauf=tauf,U=Ustp)
  mystpEEN=STP(C2,taud=taud,tauf=tauf,U=Ustp)
  #########
  #########
  #########
  spikes=SpikeMonitor(networkE)
  ### 1st step simulation: until stim on
  run(stim_on,report='text')
  ##2nd step of the simulation: stim presentation
  ## If I am not wrong, here I specify the location with this "float(NE)-0.5", being on the center
  ## pos_stim is a fraction where 0 is the first neuron (0) and 1 is the last one(360). 0.5 means stimulation in the middle (180)
  ## 0.75 means stimulating at 90 (1-0.75=0.25; 0.35*360=90). 0.25 means stimulating at 270 (1-0.25=0.75; 0.75*360=270)
  pos=arange(NE)
  networkE.Iext=stimE*exp(-0.5*(pos/float(NE)-pos_stim)**2/(epsE**2))#stimE*(1.+epsE*cos(2*pi*(pos/float(NE)-0.5)))
  pos=arange(NI)
  networkI.Iext=stimI*(1.+epsI*cos(2*pi*(pos/float(NI)-pos_stim)))
  run(stim_off-stim_on,report='text')
  ##3rd step of the simulation: delay period
  networkE.Iext=extE
  networkI.Iext=extI
  counts=SpikeCounter(networkE)
  run(runtime-stim_off,report='text')
  #last_sec=10*second - 1*second
  #rates=counts.count/(runtime-last_sec)
  return spikes.it


