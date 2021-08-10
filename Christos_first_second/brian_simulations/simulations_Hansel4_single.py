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
print('0')
one_simulation = run_simulation(IEext=0., pos_stim=0.5, paralel_simulation=False, path_='/home/david/Desktop/brian_simulations_albert/')

print('1.25')
one_simulation = run_simulation(IEext=1.25, pos_stim=0.5, paralel_simulation=False, path_='/home/david/Desktop/brian_simulations_albert/')



