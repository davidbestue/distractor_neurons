# -*- coding: utf-8 -*-
"""
@author: David Bestue
"""


from model_phantom_DB import *
numcores = multiprocessing.cpu_count() - 3
from joblib import Parallel, delayed
import multiprocessing


N=512
ch_size = 90
ch = int(360/ch_size)
reps_=3



rEs_on = []
targets = []

for iPos in np.arange(0, 360, ch_size):
    for n in range(reps_):
        on_far_X= model(totalTime=1000, targ_onset_1=0, targ_onset_2=4000, angle_target_i=iPos, presentation_period=100,
                   angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=0.05, I0I=0.5, 
                   GEE=0.068*fee,
                   GII= 0.13*fii,
                   GEI=0.13*fei,
                   GIE=0.042*fie, 
                   sigE=7., sigI=5., k_noise=0.6,           
                   kappa_E=45, 
                   kappa_I=0.3, 
                   kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
                   plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
                   phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)
        targets.append(iPos)
        rEs_on.append(on_far_X[3])










paths_save_= '/home/david/Desktop/IDIBAPS/Simulations_radial/results_simulations_radial_linear_noiser2.xlsx'

frames=[]

for idx, TIMES in enumerate(list(np.arange(0,4000, 1000) + 450 ) ): ##4000
	print(TIMES)
	Positions = list(np.arange(2,4.75,0.25))*500 ##0.25
	Times=[TIMES for i in range(len(Positions))]
	outputs= Parallel(n_jobs = numcores)(delayed(model_radial_linear)(totalTime=tim, 
	           targ_onset=100,  
	           presentation_period=350,
	           position=posx, 
	           tauE=9, tauI=4,  
	           I0E=0.1, I0I=0.5,
	           GEE=0.022, GEI=0.019, GIE=0.01 , GII=0.1, 
	           NsigE=0.95, NsigI=1.8, 
	           N=512, rint = 1, rext = 6,
	           plot_connectivity=False, 
	           plot_rate=False, save_RE=False) for posx, tim in zip(Positions, Times)) 
	#
	df = pd.DataFrame(outputs)
	df.columns=['interference', 'position', 'simul_time']
	df['delay_time']=TIMES-450
	frames.append(df)
	############


##
df_tot= pd.concat(frames)
df_tot.to_excel(paths_save_)


