from model_phantom import *
#from linares_plot import * 

n_simuls=500 # 50
numcores = multiprocessing.cpu_count() -5
print('Number cores: '+ str(numcores))



fee=1.09
fei=1.27
fie=0.96
fii=0.33



ON_1_far = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=300, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 



OFF_1_far = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=300, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-6., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 




err1_on_f = pd.DataFrame([ON_1_far[i][1] for i in range(len(ON_1_far))])
err1_on_f.columns=['err']
err1_on_f['abs_err']=abs(err1_on_f['err'])
err1_on_f['stimulation']='ON'
err1_on_f['distance']='far'
err1_on_f['order']='1st'
err1_on_f.to_excel('/home/david/Desktop/err1_on_f24.xlsx') 



err1_off_f = pd.DataFrame([OFF_1_far[i][1] for i in range(len(OFF_1_far))])
err1_off_f.columns=['err']
err1_off_f['abs_err']=abs(err1_off_f['err'])
err1_off_f['stimulation']='OFF'
err1_off_f['distance']='far'
err1_off_f['order']='1st'
err1_off_f.to_excel('/home/david/Desktop/err1_off_f24.xlsx')





ON_1_close = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=62, tauE=20, tauI=10,  n_stims=2, I0E=-1, I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 




OFF_1_close =Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=62, tauE=20, tauI=10,  n_stims=2, I0E=-6., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 




err1_on_c = pd.DataFrame([ON_1_close[i][1] for i in range(len(ON_1_close))])
err1_on_c.columns=['err']
err1_on_c['abs_err']=abs(err1_on_c['err'])
err1_on_c['stimulation']='ON'
err1_on_c['distance']='close'
err1_on_c['order']='1st'
err1_on_c.to_excel('/home/david/Desktop/err1_on_c23.xlsx')


err1_off_c = pd.DataFrame([OFF_1_close[i][1] for i in range(len(OFF_1_close))])
err1_off_c.columns=['err']
err1_off_c['abs_err']=abs(err1_off_c['err'])
err1_off_c['stimulation']='OFF'
err1_off_c['distance']='close'
err1_off_c['order']='1st'
err1_off_c.to_excel('/home/david/Desktop/err1_off_c23.xlsx')




fee=1.088
fei=1.045
fie=1.42
fii=1.01



ON_2_far = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-1., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 



OFF_2_far = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-6., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 



err2_on_f = pd.DataFrame([ON_2_far[i][2] for i in range(len(ON_2_far))])
err2_on_f.columns=['err']
err2_on_f['abs_err']=abs(err2_on_f['err'])
err2_on_f['stimulation']='ON'
err2_on_f['distance']='far'
err2_on_f['order']='2nd'
err2_on_f.to_excel('/home/david/Desktop/err2_on_f23.xlsx')




err2_off_f = pd.DataFrame([OFF_2_far[i][2] for i in range(len(OFF_2_far))])
err2_off_f.columns=['err']
err2_off_f['abs_err']=abs(err2_off_f['err'])
err2_off_f['stimulation']='OFF'
err2_off_f['distance']='far'
err2_off_f['order']='2nd'
err2_off_f.to_excel('/home/david/Desktop/err2_off_f23.xlsx')



ON_2_close = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=62, tauE=20, tauI=10,  n_stims=2, I0E=-1., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 




OFF_2_close = Parallel(n_jobs = numcores)(delayed(model)(totalTime=2000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=62, tauE=20, tauI=10,  n_stims=2, I0E=-6., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fii,
           GEI=0.13*fei,
           GIE=0.042*fie, 
           sigE=3.25, sigI=2.,            
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=12.5, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=3.9, phantom_onset=50000, phnatom_duration=100)  for n in range(n_simuls)) 





err2_on_c = pd.DataFrame([ON_2_close[i][2] for i in range(len(ON_2_close))])
err2_on_c.columns=['err']
err2_on_c['abs_err']=abs(err2_on_c['err'])
err2_on_c['stimulation']='ON'
err2_on_c['distance']='close'
err2_on_c['order']='2nd'
err2_on_c.to_excel('/home/david/Desktop/err2_on_c23.xlsx')


err2_off_c = pd.DataFrame([OFF_2_close[i][2] for i in range(len(OFF_2_close))])
err2_off_c.columns=['err']
err2_off_c['abs_err']=abs(err2_off_c['err'])
err2_off_c['stimulation']='OFF'
err2_off_c['distance']='close'
err2_off_c['order']='2nd'
err2_off_c.to_excel('/home/david/Desktop/err2_off_c23.xlsx')





