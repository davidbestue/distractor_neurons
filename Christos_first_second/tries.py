from model import *

numcores = multiprocessing.cpu_count() -1 
print('Number cores: '+ str(numcores))
n_simuls=100

## ### Remember 1st
fee=1. 
fei=1. 
fie=1. 
fii=1.


### Close: off
results_1st_close_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
           angle_separation=51, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #OFF
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 



OFF= pd.DataFrame( [results_1st_close_off[i][1] for i in range(len(results_1st_close_off))])
OFF['stimul']='OFF' 
OFF['position']='close' 


### Close: on
results_1st_close_on = Parallel(n_jobs = numcores)(delayed(model)(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
           angle_separation=51, tauE=20, tauI=10,  n_stims=2, I0E=0.18, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #OFF
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 



ON= pd.DataFrame( [results_1st_close_on[i][1] for i in range(len(results_1st_close_on))])
ON['stimul']='ON' 
ON['position']='close' 

##
df = pd.concat([OFF, ON], ignore_index=True)