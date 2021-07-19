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
# one_simulation =  simulation(loadconnections=False, saveconnections=True, save_name='connections_sp_20000', N=20000)
# io.savemat('/home/david/Desktop/brian_simulations/single_simulation_20000',{'spktm': one_simulation})


## One simulation
one_simulation =  simulation(name_conections='connections_sp_20000.npz', N=20000)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation2_20000',{'spktm': one_simulation})



