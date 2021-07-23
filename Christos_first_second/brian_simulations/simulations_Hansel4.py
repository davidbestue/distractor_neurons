#import brian_no_units


from model_Hansel4 import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3



## One simulation
#one_simulation = run_simulation(IEext=0.2, pos_stim=0.75, save_file=True)
##one_simulation = run_simulation(IEext=0.2, pos_stim=0.25, save_file=False)

#######################################


### Multiple simulations in paralel

extEs = [0, 0.2, 0.5, 1, 2]
n_ext = len(extEs)
positions = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
n_pos = len(positions)

number_ = 100

extEs_p = extEs*n_pos*number_
positions_p = positions*n_ext*number_ 

results = Parallel(n_jobs = numcores)(delayed(run_simulation)(IEext=extE, pos_stim=pos, save_file=False)  for extE, pos in zip(extEs_p, positions_p))    

# io.savemat('/home/david/Desktop/brian_simulations/results_simulations_1000',{'extEs':extEs, 'spktm': results})





######################################



# start_time = time.time()
# numcores   = mp.cpu_count()
# defaultclock.reinit()
# defaultclock.dt = 0.1*ms


# IEext=0.
# pos_stim=0.5
# save_file=False
# i=1


# global par
# time_s = int(str(time.time()).split('.')[0])
# save_name    = "simulation_%i_%i_%s" %(os.getpid(), time_s, socket.gethostname())
# print save_name

# loadconnections = True #False
# saveconnections = False #True

# ### nstims=50
# stim_on=2000*ms
# stim_off=3000*ms
# stimE=0.1*2.4*mV # stimulus input
# epsE=60
# stimI=0*mV
# epsI=0
# runtime=7*second

# N=20000 #10000 # total number of neurons
# K=500 #250 # total number of inputs
# tauE=20*ms 
# tauI=10*ms 

# taua = 3*ms # AMPA synapse decay
# taun = 50*ms # NMDA synapse decay
# taug = 4*ms # GABA synapse decay
# Ustp = 0.04 #0.03
# taud = 200*ms
# tauf = 450*ms
# R0 = 20.0/second #28.0/second #peak imbl firing rate 28Hz, mean imbl firing rate 14 Hz
# u_ss = Ustp*(1.+R0*tauf)/(1.+Ustp*R0*tauf)
# u_ss0= 0.03*(1.+R0*tauf)/(1.+0.03*R0*tauf)
# synfact_ss = u_ss/(1.+taud*R0*u_ss)
# synfact_ss0 = u_ss0/(1.+taud*R0*u_ss0)

# Vt  = 20*mV          # spike threshold
# Vr  = -3.33*mV          # reset value
# refE= 0*ms                # refractory periods
# refI= 0*ms                # refractory periods

# gEEA=533.3*mV*ms  
# gEEN=0.92*533.3*mV*ms  
# gEIA=67.2*mV*ms  
# gEIN=1*7.4*mV*ms  #imbl: factor 0.5
# gIE=-138.6*mV*ms
# gII=-90.6*mV*ms
# sigmaEE=30 #60 # E-to-E footprint in degrees
# sigmaEI=35 #70 # E-to-E footprint in degrees
# sigmaIE=30 #60 # E-to-E footprint in degrees
# sigmaII=30 #60 # E-to-E footprint in degrees

# #these are intermediate calculations needed for the equations below

# NE=int(ceil(0.8*N)) # number of excitatory neurons
# NI=int(floor(0.2*N)) # number of inhibitory neurons

# KE=int(ceil(0.8*K)) # number of excitatory inputs
# KI=int(floor(0.2*K)) # number of inhibitory inputs

# sigEE=sigmaEE/360.0
# sigEI=sigmaEI/360.0
# sigIE=sigmaIE/360.0
# sigII=sigmaII/360.0
# epsE = epsE/360.0

# gEEA=gEEA/sqrt(KE)/taua*synfact_ss0/synfact_ss #
# gEEN=gEEN/sqrt(KE)/taun*synfact_ss0/synfact_ss #
# gEIA=gEIA/sqrt(KE)/taua
# gEIN=gEIN/sqrt(KE)/taun
# gIE=gIE/sqrt(KI)/taug
# gII=gII/sqrt(KI)/taug

# stimE=stimE*sqrt(KE)
# stimI=stimI*sqrt(KE)

# #ring model connectivity
# def conn(k,sig,C):
#     if abs(k)<=0.5:
#         value=C*exp(-0.5*(k/sig)**2)/sqrt(2*pi)/sig
#     else:
#         value=C*exp(-0.5*((abs(k)-1)/sig)**2)/sqrt(2*pi)/sig
#     return value

# #equations for each neuron. xpre and satura are shadow variables that track the saturated state
# # of a neuron's outgoing NMDA synapses, and this is used as a "synaptic depression" effect on x
# # (see modulation='satura' below) so as to mimic the NMDA saturation that is modeled in 
# # Compte et al., 2000
# eqs = '''
# dV/dt = (gea+gen+gi-V+Ix+Iext)/tau : mV
# dgea/dt = -gea/taua          : mV
# dgen/dt = -gen/taun   : mV
# dgi/dt = -gi/taug             : mV
# tau : second
# Ix : mV
# Iext : mV
# '''

# networkE=NeuronGroup(NE,model=eqs,threshold=Vt,reset=Vr, refractory=refE)
# networkI=NeuronGroup(NI,model=eqs,threshold=Vt,reset=Vr, refractory=refI)
# networkE.tau=tauE
# networkI.tau=tauI
# networkE.Ix=1.66*sqrt(KE)*mV
# networkI.Ix=1.85*0.83*sqrt(KE)*mV
# networkE.V =  Vr + rand(NE)*(Vt - Vr) #Vt-2.0*mV + rand(NE) * 2.0*mV
# networkI.V =  Vr + rand(NI)*(Vt - Vr) #Vt-2.0*mV + rand(NI) * 2.0*mV

# networkE.Iext=IEext*mV
# networkI.Iext=0*mV

# if loadconnections:
#     loader = np.load('connections_sp.npz', allow_pickle=True)
#     WW = sparse.csr_matrix((loader['CEEd'], loader['CEEi'], loader['CEEp']), shape=loader['CEEs'])
#     WW[WW != 0] = 1
#     C1=Connection(networkE, networkE, 'gea', weight=gEEA*WW)
#     C2=Connection(networkE, networkE, 'gen', weight=gEEN*WW)
#     WW = sparse.csr_matrix((loader['CEId'], loader['CEIi'], loader['CEIp']), shape=loader['CEIs'])
#     WW[WW != 0] = 1
#     C3=Connection(networkE, networkI, 'gea', weight=gEIA*WW)
#     C4=Connection(networkE, networkI, 'gen', weight=gEIN*WW)
#     WW = sparse.csr_matrix((loader['CIEd'], loader['CIEi'], loader['CIEp']), shape=loader['CIEs'])
#     WW[WW != 0] = 1
#     C5=Connection(networkI, networkE, 'gi', weight=gIE*WW)
#     WW = sparse.csr_matrix((loader['CIId'], loader['CIIi'], loader['CIIp']), shape=loader['CIIs'])
#     WW[WW != 0] = 1
#     C6=Connection(networkI, networkI, 'gi', weight=gII*WW)
# else:
#     fE = float(KE)/float(NE)
#     fI = float(KI)/float(NI)
#     C1=Connection(networkE, networkE, 'gea', weight=gEEA, sparseness=lambda ix,j: conn(float(ix)/NE-float(j)/NE,sigEE,fE))
#     C2=Connection(networkE, networkE, 'gen', weight=gEEN/gEEA*C1.W)
#     C3=Connection(networkE, networkI, 'gea', weight=gEIA, sparseness=lambda ix,j: conn(float(ix)/NE-float(j)/NI,sigEI,fE))
#     C4=Connection(networkE, networkI, 'gen', weight=gEIN/gEIA*C3.W)
#     C5=Connection(networkI, networkE, 'gi', weight=gIE, sparseness=lambda ix,j: conn(float(ix)/NI-float(j)/NE,sigIE,fI))
#     C6=Connection(networkI, networkI, 'gi', weight=gII, sparseness=lambda ix,j: conn(float(ix)/NI-float(j)/NI,sigII,fI))

# if saveconnections:
#     CEE = sparse.csr_matrix(C1.W)
#     CEI = sparse.csr_matrix(C3.W)
#     CIE = sparse.csr_matrix(C5.W)
#     CII = sparse.csr_matrix(C6.W)
#     np.savez_compressed('connections_sp', CEEd=CEE.data, CEEi=CEE.indices, CEEp=CEE.indptr, CEEs=CEE.shape,
#                                          CEId=CEI.data, CEIi=CEI.indices, CEIp=CEI.indptr, CEIs=CEI.shape,
#                                          CIEd=CIE.data, CIEi=CIE.indices, CIEp=CIE.indptr, CIEs=CIE.shape,
#                                          CIId=CII.data, CIIi=CII.indices, CIIp=CII.indptr, CIIs=CII.shape)


# mystpEEA=STP(C1,taud=taud,tauf=tauf,U=Ustp)
# mystpEEN=STP(C2,taud=taud,tauf=tauf,U=Ustp)

# spikes=SpikeMonitor(networkE)

# run(stim_on,report='text')

# #  
# ## pos_stim is a fraction where 0 is the first neuron (0) and 1 is the last one(360). 0.5 means stimulation in the middle (180)
# ## 0.75 means stimulating at 90 (1-0.75=0.25; 0.35*360=90). 0.25 means stimulating at 270 (1-0.25=0.75; 0.75*360=270)
# pos=arange(NE)
# networkE.Iext=stimE*exp(-0.5*(pos/float(NE)-pos_stim)**2/(epsE**2))#stimE*(1.+epsE*cos(2*pi*(pos/float(NE)-0.5)))
# pos=arange(NI)
# networkI.Iext=stimI*(1.+epsI*cos(2*pi*(pos/float(NI)-pos_stim)))
# #
# run(stim_off-stim_on,report='text')

# networkE.Iext=IEext*mV
# networkI.Iext=0*mV

# counts=SpikeCounter(networkE)

# run(runtime-stim_off,report='text')

# rates=counts.count/(runtime-stim_off)

# i,t         = spikes.it

# popdectm, nwins = readout(i, t, runtime, NE)
# dectm = [np.degrees(popdectm[i]) for i in range(len(popdectm))] 
# errdec = [np.degrees(popdectm[i])-360*(1-pos_stim) for i in range(len(popdectm))] 

