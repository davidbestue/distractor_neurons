#import brian_no_units


from model import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3



## One simulation
one_simulation =  simulation(name_conections='connections_sp_1000.npz', N=1000, extE=0.5, pos_stim=0.75)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation_90deg',{'spktm': one_simulation})




