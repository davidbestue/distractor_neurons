#import brian_no_units


from model_Hansel4 import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3



## One simulation

run_simulation(IEext=0, pos_stim=0.5, paralel_simulation=False, path_='/home/david/Desktop/brian_simulations_albert/')
#one_simulation = run_simulation(IEext=0.2, pos_stim=0.75, save_file=True)
##one_simulation = run_simulation(IEext=0.2, pos_stim=0.25, save_file=False)

#######################################


### Multiple simulations in paralel

# extEs = [0, 0.5, 1]
# #n_ext = len(extEs)
# positions = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# #n_pos = len(positions)

# number_ = 5


# for extE in extEs:
# 	for pos in positions:
# 		for n in range(number_):
# 			results = run_simulation(IEext=extE, pos_stim=pos, paralel_simulation=True, path_='/home/david/Desktop/brian_simulations_albert/simulations2')





