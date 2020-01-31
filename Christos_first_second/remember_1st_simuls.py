
from model import *

## Regime for remember 1st
fee=1. 
fei=1. 
fie=1. 
fii=1.


### Example
# an = model(totalTime=1500, targ_onset_1=100, targ_onset_2=800, angle_target_i=90, presentation_period=100,
#            angle_separation=100, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.12*fie,
#            GIE=0.042*fii, 
#            sigE=10., sigI=5., 
#            kappa_E=45, 
#            kappa_I=0.5, 
#            kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False) 


n_simuls=5
# numcores = multiprocessing.cpu_count() 
# print('Number cores: '+ str(numcores))
# results_1st_close_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
#            angle_separation=49, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=10., sigI=5., 
#            kappa_E=45, 
#            kappa_I=0.5, #OFF
#            kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 


### Close: off
res_off=[]

for i in range(n_simuls):
    an = model(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
               angle_separation=49, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
               GEE=0.068*fee,
               GII= 0.13*fei,
               GEI=0.13*fie,
               GIE=0.042*fii, 
               sigE=10., sigI=5., 
               kappa_E=45, 
               kappa_I=0.5, #OFF
               kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
               plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False) 

    print(an[1], an[2])
    res_off.append(an[1])


###
print('abs error OFF, close: ' + str(round(np.mean([abs(res_off[i]) for i in range(len(res_off))]),2) ) )

OFF = pd.DataFrame(res_off)
OFF['stimul']='ON' 


### Close: on
res_on=[]

for i in range(n_simuls):
    an = model(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
               angle_separation=49, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
               GEE=0.068*fee,
               GII= 0.13*fei,
               GEI=0.13*fie,
               GIE=0.042*fii, 
               sigE=10., sigI=5., 
               kappa_E=45, 
               kappa_I=0.35, #OFF
               kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
               plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False) 

    print(an[1], an[2])
    res_on.append(an[1])


###
print('abs error ON, close: ' + str(round(np.mean([abs(res_on[i]) for i in range(len(res_on))]),2) ) )



ON = pd.DataFrame(res_on)
ON['stimul']='ON' 