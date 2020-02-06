
from model import *

numcores = multiprocessing.cpu_count() - 1
print('Number cores: '+ str(numcores)) 
n_simuls=400

## Remember 1st

fee=1. 
fei=1. 
fie=1. 
fii=1.


### Far: off
results_1st_far_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
           angle_separation=135, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #OFF
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

#np.mean([abs(results_1st_far_off[i][1]) for i in range(len(results_1st_far_off))]) 
OFF_f= pd.DataFrame( [results_1st_far_off[i][2] for i in range(len(results_1st_far_off))])
OFF_f['stimul']='OFF' 
OFF_f['position']='far'  

### Far: on
results_1st_far_on = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
           angle_separation=135, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.35, #ON
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

#np.mean([abs(results_1st_far_on[i][1]) for i in range(len(results_1st_far_on))]) 
ON_f= pd.DataFrame( [results_1st_far_on[i][2] for i in range(len(results_1st_far_on))])
ON_f['stimul']='ON' 
ON_f['position']='far' 


df_f = pd.concat([OFF_f, ON_f], ignore_index=True)
###df_f.to_excel('/home/david/Desktop/remembers_second_far.xlsx')
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

#df_ =pd.concat([df_c, df_f], ignore_index=True)
df_f.to_excel('/home/david/Desktop/135_remember_first.xlsx')

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################



#### Borrar esto despues!!!!

fee=0.94
fei=0.92
fie=1.14
fii=1.08


### Far: off
results_2nd_far_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
           angle_separation=135, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #OFF
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

#np.mean([abs(results_2nd_far_off[i][1]) for i in range(len(results_2nd_far_off))]) 
OFF_f= pd.DataFrame( [results_2nd_far_off[i][2] for i in range(len(results_2nd_far_off))])
OFF_f['stimul']='OFF' 
OFF_f['position']='far'  

### Far: on
results_2nd_far_on = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
           angle_separation=135, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.35, #ON
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

#np.mean([abs(results_2nd_far_on[i][1]) for i in range(len(results_2nd_far_on))]) 
ON_f= pd.DataFrame( [results_2nd_far_on[i][2] for i in range(len(results_2nd_far_on))])
ON_f['stimul']='ON' 
ON_f['position']='far' 


df_f = pd.concat([OFF_f, ON_f], ignore_index=True)
###df_f.to_excel('/home/david/Desktop/remembers_second_far.xlsx')
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

#df_ =pd.concat([df_c, df_f], ignore_index=True)
df_f.to_excel('/home/david/Desktop/135_remember_second.xlsx')



