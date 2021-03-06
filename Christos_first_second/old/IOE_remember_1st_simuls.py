
from model import *

numcores = multiprocessing.cpu_count() -1 
print('Number cores: '+ str(numcores))
n_simuls=400

## ### Remember 1st
fee=1. 
fei=1. 
fie=1. 
fii=1.

# an = model(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
#            angle_separation=51, tauE=20, tauI=10,  n_stims=2, I0E=0.18, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=10., sigI=5., 
#            kappa_E=45, 
#            kappa_I=0.5, 
#            kappa_stim=40., N=512, stim_strengthE=9.2, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False) 

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

#ON=pd.concat([ON, ON2], ignore_index=True)
##
df = pd.concat([OFF, ON], ignore_index=True)
#df.to_excel('/home/david/Desktop/remembers_first_close_51_del.xlsx')


### Far: off

results_1st_far_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
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

OFF_f= pd.DataFrame( [results_1st_far_off[i][1] for i in range(len(results_1st_far_off))])
OFF_f['stimul']='OFF' 
OFF_f['position']='far'  

### Far: on
results_1st_far_on = Parallel(n_jobs = numcores)(delayed(model)(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
           angle_separation=135, tauE=20, tauI=10,  n_stims=2, I0E=0.18, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #ON
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

ON_f= pd.DataFrame( [results_1st_far_on[i][1] for i in range(len(results_1st_far_on))])
ON_f['stimul']='ON' 
ON_f['position']='far' 


df_f = pd.concat([OFF_f, ON_f], ignore_index=True)
###df_f.to_excel('/home/david/Desktop/remembers_first_far.xlsx')

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

df_ =pd.concat([df, df_f], ignore_index=True)
df_.to_excel('/home/david/Desktop/remembers_first_I0E_del.xlsx')

###df_.to_excel('remembers_first_I0E_del.xlsx')

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

## #serie
## res_off=[]
## for i in range(n_simuls):
##     an = model(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
##                angle_separation=49, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
##                GEE=0.068*fee,
##                GII= 0.13*fei,
##                GEI=0.13*fie,
##                GIE=0.042*fii, 
##                sigE=10., sigI=5., 
##                kappa_E=45, 
##                kappa_I=0.5, #OFF
##                kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
##                plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False) 

##     print(an[1], an[2])
##     res_off.append(an[1])
###



## ### Remember 2nd
fee=0.94
fei=0.92
fie=1.14
fii=1.08

# an = model(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
#            angle_separation=51, tauE=20, tauI=10,  n_stims=2, I0E=0.18, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=10., sigI=5., 
#            kappa_E=45, 
#            kappa_I=0.5, 
#            kappa_stim=40., N=512, stim_strengthE=9.2, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False) 

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



OFF= pd.DataFrame( [results_1st_close_off[i][2] for i in range(len(results_1st_close_off))])
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
           kappa_I=0.5, #ON
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

#
ON= pd.DataFrame( [results_1st_close_on[i][2] for i in range(len(results_1st_close_on))])
ON['stimul']='ON' 
ON['position']='close' 

#ON=pd.concat([ON, ON2], ignore_index=True)
##
df = pd.concat([OFF, ON], ignore_index=True)
#df.to_excel('/home/david/Desktop/remembers_first_close_51_del.xlsx')


### Far: off

results_1st_far_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
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

OFF_f= pd.DataFrame( [results_1st_far_off[i][2] for i in range(len(results_1st_far_off))])
OFF_f['stimul']='OFF' 
OFF_f['position']='far'  

### Far: on
results_1st_far_on = Parallel(n_jobs = numcores)(delayed(model)(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
           angle_separation=135, tauE=20, tauI=10,  n_stims=2, I0E=0.18, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #ON
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

ON_f= pd.DataFrame( [results_1st_far_on[i][2] for i in range(len(results_1st_far_on))])
ON_f['stimul']='ON' 
ON_f['position']='far' 


df_f = pd.concat([OFF_f, ON_f], ignore_index=True)
###df_f.to_excel('/home/david/Desktop/remembers_first_far.xlsx')

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

df_ =pd.concat([df, df_f], ignore_index=True)
df_.to_excel('/home/david/Desktop/remembers_second_I0E_del.xlsx')



