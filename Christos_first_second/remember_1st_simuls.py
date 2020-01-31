
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


n_simuls=12
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



res=[]

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
    res.append(an[1])


###

print(np.mean(res))
