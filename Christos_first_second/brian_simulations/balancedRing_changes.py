#import brian_no_units

from brian import *
import scipy.io as io
import cPickle
import sys
from StringIO import StringIO
from scipy import sparsenes
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

defaultclock.reinit()
defaultclock.dt = 0.1*ms

loadconnections = False #True
saveconnections =  True #False

stim_on=2000*ms
stim_off=3000*ms
stimE=0.1*2.4*mV # stimulus input
epsE=30
stimI=0*mV
epsI=0
runtime=10*second

N=1000 ##30000 # total number of neurons
K=1500 #250 # total number of inputs
tauE=20*ms 
tauI=10*ms 

taua = 3*ms # AMPA synapse decay
taun = 50*ms # NMDA synapse decay
taug = 4*ms # GABA synapse decay
Ustp = 0.03 # 0.03
taud = 200*ms
tauf = 450*ms
R0 = 28.0/second
u_ss = Ustp*(1.+R0*tauf)/(1.+Ustp*R0*tauf)
u_ss0= 0.03*(1.+R0*tauf)/(1.+0.03*R0*tauf)
synfact_ss = u_ss/(1.+taud*R0*u_ss)
synfact_ss0 = u_ss0/(1.+taud*R0*u_ss0)
ntaud=(u_ss-synfact_ss0)/(synfact_ss0*R0*u_ss) #to avoid scaling

Vt  = 20*mV          # spike threshold
Vr  = -3.33*mV          # reset value
refE= 0*ms                # refractory periods
refI= 0*ms                # refractory periods

gEEA=533.3*mV*ms  
gEEN=0.95*533.3*mV*ms  
gEIA=67.2*mV*ms  
gEIN=1*7.4*mV*ms
gIE=-138.6*mV*ms
gII=-90.6*mV*ms
sigmaEE=30 #60 # E-to-E footprint in degrees
sigmaEI=35 #70 # E-to-E footprint in degrees
sigmaIE=30 # 60 # E-to-E footprint in degrees
sigmaII=30 #60 # E-to-E footprint in degrees

#these are intermediate calculations needed for the equations below

NE=int(ceil(0.8*N)) # number of excitatory neurons
NI=int(floor(0.2*N)) # number of inhibitory neurons

KE=int(ceil(0.8*K)) # number of excitatory inputs
KI=int(floor(0.2*K)) # number of inhibitory inputs

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

#ring model connectivity
def conn(k,sig,C):
  if abs(k)<=0.5:
    value=C*exp(-0.5*(k/sig)**2)/sqrt(2*pi)/sig
  else:
    value=C*exp(-0.5*((abs(k)-1)/sig)**2)/sqrt(2*pi)/sig
  return value

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

networkE=NeuronGroup(NE,model=eqs,threshold=Vt,reset=Vr, refractory=refE)
networkI=NeuronGroup(NI,model=eqs,threshold=Vt,reset=Vr, refractory=refI)
networkE.tau=tauE
networkI.tau=tauI
networkE.Ix=1.66*sqrt(KE)*mV
networkI.Ix=1.85*0.83*sqrt(KE)*mV
networkE.V = Vr + rand(NE)*(Vt - Vr) #Vt-2.0*mV + rand(NE) * 2.0*mV
networkI.V = Vr + rand(NI)*(Vt - Vr) #Vt-2.0*mV + rand(NI) * 2.0*mV

networkE.Iext=0*mV
networkI.Iext=0*mV


fE = float(KE)/float(NE)
fI = float(KI)/float(NI)

C1=Connection(networkE, networkE, 'gea', weight=gEEA, sparseness=lambda i,j: conn(float(i)/NE-float(j)/NE,sigEE,fE))
C2=Connection(networkE, networkE, 'gen', weight=gEEN/gEEA*C1.W)
C3=Connection(networkE, networkI, 'gea', weight=gEIA, sparseness=lambda i,j: conn(float(i)/NE-float(j)/NI,sigEI,fE))
C4=Connection(networkE, networkI, 'gen', weight=gEIN/gEIA*C3.W)
C5=Connection(networkI, networkE, 'gi', weight=gIE, sparseness=lambda i,j: conn(float(i)/NI-float(j)/NE,sigIE,fI))
C6=Connection(networkI, networkI, 'gi', weight=gII, sparseness=lambda i,j: conn(float(i)/NI-float(j)/NI,sigII,fI))




  CEE = csr_matrix(C1.W)
  CEI = csr_matrix(C3.W)
  CIE = csr_matrix(C5.W)
  CII = csr_matrix(C6.W)
  np.savez_compressed('connections_sp', CEEd=CEE.data, CEEi=CEE.indices, CEEp=CEE.indptr, CEEs=CEE.shape, 
            CEId=CEI.data, CEIi=CEI.indices, CEIp=CEI.indptr, CEIs=CEI.shape, 
            CIEd=CIE.data, CIEi=CIE.indices, CIEp=CIE.indptr, CIEs=CIE.shape, 
            CIId=CII.data, CIIi=CII.indices, CIIp=CII.indptr, CIIs=CII.shape)
#  f_out = open('connections','w') #
#  cPickle.dump(C1.W,f_out)
#  cPickle.dump(C2.W,f_out)
#  cPickle.dump(C3.W,f_out)
#  cPickle.dump(C4.W,f_out)
#  cPickle.dump(C5.W,f_out)
#  cPickle.dump(C6.W,f_out)
#  f_out.close()

mystpEEA=STP(C1,taud=taud,tauf=tauf,U=Ustp)
mystpEEN=STP(C2,taud=taud,tauf=tauf,U=Ustp)

spikes=SpikeMonitor(networkE)

run(stim_on,report='text')

pos=arange(NE)
networkE.Iext=stimE*exp(-0.5*(pos/float(NE)-0.5)**2/(epsE**2))#stimE*(1.+epsE*cos(2*pi*(pos/float(NE)-0.5)))
pos=arange(NI)
networkI.Iext=stimI*(1.+epsI*cos(2*pi*(pos/float(NI)-0.5)))

run(stim_off-stim_on,report='text')

networkE.Iext=0*mV
networkI.Iext=0*mV

counts=SpikeCounter(networkE)

run(runtime-stim_off,report='text')

rates=counts.count/(runtime-stim_off)

#io.savemat('results_balancedRing',{'rate':rates, 'spktm': spikes.it})

plt.plot(spikes.it[1],spikes.it[0],'k.',markersize=2)
plt.show()


