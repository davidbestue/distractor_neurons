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



def simulation(loadconnections=True, name_conections='connections_sp_30000.npz', 
  saveconnections=False, save_name='connections_sp_30000',
  dt_clock=0.1, stimon=2000, stimoff=3000,  
  epsE=30, epsI=0, timesimulation=10, N=30000, prop_e=0.8, prop_i=0.2, K=1500, 
  tE = 20, tI=10, ta=3, tn=50, tg=4, td=200, tf=450, 
  Ustp = 0.03, ro=28.0, vt=20, vr=-3.33, fe=0., fi=0., 
  eea=533.3, een=0.95, eia=67.2, ein=7.4, ie=-138.6, ii=-90.6, 
  sigmaEE=30, sigmaEI=35, sigmaIE=30, sigmaII=30, 
  extE=0., extI=0., stim_str = 0.24, pos_stim=0.5):
    #############
    defaultclock.reinit()
    defaultclock.dt = dt_clock*ms
    stim_on=stimon*ms
    stim_off=stimoff*ms
    stimE=stim_str*mV # stimulus input
    stimI=0*mV
    runtime=timesimulation*second
    ### Taus for the spiking rate model of persistent activity
    tauE=tE*ms 
    tauI=tI*ms 
    taua = ta*ms # AMPA synapse decay
    taun = tn*ms # NMDA synapse decay
    taug = tg*ms # GABA synapse decay
    ### STP equations (http://www.scholarpedia.org/article/Short-term_synaptic_plasticity)
    taud = td*ms # synaptic depression tau
    tauf = tf*ms # synaptic facilitation tau
    R0 = ro/second
    u_ss = Ustp*(1.+R0*tauf)/(1.+Ustp*R0*tauf)
    u_ss0= 0.03*(1.+R0*tauf)/(1.+0.03*R0*tauf)
    synfact_ss = u_ss/(1.+taud*R0*u_ss)
    synfact_ss0 = u_ss0/(1.+taud*R0*u_ss0)
    ntaud=(u_ss-synfact_ss0)/(synfact_ss0*R0*u_ss) #to avoid scaling
    Vt  = vt*mV          # spike threshold
    Vr  = vr*mV          # reset value
    refE= fe*ms          # refractory periods
    refI= fi*ms          # refractory periods
    ### Conductances
    gEEA=eea*mV*ms  
    gEEN=een*eea*mV*ms  
    gEIA=eia*mV*ms  
    gEIN=ein*mV*ms
    gIE=ie*mV*ms
    gII=ii*mV*ms
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
    networkE.Iext=extE*mV #### External input (IE0 en el rate model)
    networkI.Iext=extI*mV
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
    networkE.Iext=extE*mV
    networkI.Iext=extI*mV
    counts=SpikeCounter(networkE)
    run(runtime-stim_off,report='text')
    #last_sec=10*second - 1*second
    #rates=counts.count/(runtime-last_sec)
    return spikes.it



# #### rates=counts.count/(runtime-stim_off)
# #### Save the files
# io.savemat('/home/david/Desktop/brian_simulations/results_simulation_30000',{'rate':rates, 'spktm': spikes.it})
# ##spiketime in the dictionary format
# dict_spiketimes = spikes.spiketimes
# pickle.dump( dict_spiketimes, open( "/home/david/Desktop/brian_simulations/dict_spiketimes_30000.pkl", "wb" ) )



