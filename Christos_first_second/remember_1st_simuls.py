
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


n_simuls=200




### Close: off
## Paralel
numcores = multiprocessing.cpu_count() 
print('Number cores: '+ str(numcores))
results_1st_close_off = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
           angle_separation=49, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.5, #OFF
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

np.mean([abs(results_1st_close_off[i][1]) for i in range(len(results_1st_close_off))]) ###¿¿¿16??? 15.6... okey....


#serie
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
print('abs error OFF, close: ' + str(round(np.mean([abs(res_off[i]) for i in range(len(res_off))]),2) ) ) ## 16

OFF = pd.DataFrame(res_off)
OFF['stimul']='ON' 
OFF['position']='close' 

OFF2= pd.DataFrame( [abs(results_1st_close_off[i][1]) for i in range(len(results_1st_close_off))])
OFF2['stimul']='OFF' 
OFF2['position']='close' 

OFF=pd.concat([OFF, OFF2], ignore_index=True)



### Close: on

## Paralel
numcores = multiprocessing.cpu_count() 
print('Number cores: '+ str(numcores))
results_1st_close_on = Parallel(n_jobs = numcores)(delayed(model)(totalTime=3000, targ_onset_1=50, targ_onset_2=500, angle_target_i=90, presentation_period=100,
           angle_separation=49, tauE=20, tauI=10,  n_stims=2, I0E=0.1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=10., sigI=5., 
           kappa_E=45, 
           kappa_I=0.35, #ON
           kappa_stim=40., N=512, stim_strengthE=9.20, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False)  for n in range(n_simuls)) 

np.mean([abs(results_1st_close_on[i][1]) for i in range(len(results_1st_close_on))]) ###¿¿¿18.31??? 17.8... ok...

## Serie
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
print('abs error ON, close: ' + str(round(np.mean([abs(res_on[i]) for i in range(len(res_on))]),2) ) ) ##18.31



ON = pd.DataFrame(res_on)
ON['stimul']='ON' 
ON['position']='close' 

ON2= pd.DataFrame( [abs(results_1st_close_on[i][1]) for i in range(len(results_1st_close_on))])
ON2['stimul']='ON' 
ON2['position']='close' 

ON=pd.concat([ON, ON2], ignore_index=True)


### Far: off
res_off_f=[]

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
    res_off_f.append(an[1])


###
print('abs error OFF, far: ' + str(round(np.mean([abs(res_off_f[i]) for i in range(len(res_off_f))]),2) ) )

OFF_f = pd.DataFrame(res_off_f)
OFF_f['stimul']='ON' 
OFF_f['position']='far' 


### Far: on
res_on_f=[]

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
    res_on_f.append(an[1])


###
print('abs error ON, far: ' + str(round(np.mean([abs(res_on_f[i]) for i in range(len(res_on_f))]),2) ) )



ON_f = pd.DataFrame(res_on_f)
ON_f['stimul']='ON' 
ON_f['position']='far' 



df = pd.concat([OFF, ON], ignore_index=True)
#df.to_excel('/home/david/Desktop/remembers_first_close.xlsx')

