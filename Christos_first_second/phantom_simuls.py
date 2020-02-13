from model_phantom import *
#from linares_plot import * 

n_simuls=5
numcores = multiprocessing.cpu_count() -1 
print('Number cores: '+ str(numcores))

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




fee=1
fei=1
fie=1
fii=1

OFF_1_far = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-2., I0I=0.5, 
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




err1_on_f = pd.DataFrame([ON_1_far[i][1] for i in range(len(ON_1_far))])
err1_on_f.columns=['err']
err1_on_f['abs_err']=abs(err1_on_f['err'])
err1_on_f['stimulation']='ON'
err1_on_f['distance']='far'
err1_on_f['order']='1st'
err1_on_f.to_excel('/home/david/Desktop/err1_on_f.xlsx')
err1_on_f=err1_on_f.loc[err1_on_f['abs_err']<25]


err1_off_f = pd.DataFrame([OFF_1_far[i][1] for i in range(len(OFF_1_far))])
err1_off_f.columns=['err']
err1_off_f['abs_err']=abs(err1_off_f['err'])
err1_off_f['stimulation']='OFF'
err1_off_f['distance']='far'
err1_off_f['order']='1st'
err1_off_f.to_excel('/home/david/Desktop/err1_off_f.xlsx')
err1_off_f=err1_off_f.loc[err1_off_f['abs_err']<25]

