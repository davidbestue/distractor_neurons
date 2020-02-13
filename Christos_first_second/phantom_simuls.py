from model_phantom import *
#from linares_plot import * 

n_simuls=5
numcores = multiprocessing.cpu_count() -1 
#print('Number cores: '+ str(numcores))

fee=1
fei=1
fie=1
fii=1

ON_1_far = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=0.08, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=7., sigI=5.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=2.5, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 



