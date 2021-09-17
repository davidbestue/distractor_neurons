#import brian_no_units

from model_Hansel5 import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3


## One simulation
#one_simulation = run_simulation(IEext=0.5, pos_stim=0.25, gee=480.0)

#######################################
# Multiple simulations in paralel
extEs = [0, 0.5]
n_ext = len(extEs)

n_=25
positions = [np.round(np.linspace(0,1,n_)[x],2) for x in range(n_)]

number_ = 5

for extE in extEs:
	for pos in positions:
		for n in range(number_):
			results = run_simulation(IEext=extE, pos_stim=pos, gee=480.0, paralel_simulation=False, path_='/home/david/Desktop/brian_simulations_albert/simulations8/')





##
#