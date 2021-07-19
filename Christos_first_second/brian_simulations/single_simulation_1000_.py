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
one_simulation =  simulation(loadconnections=False, saveconnections=True, save_name='connections_sp_1000_', N=1000, K=25)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation_1000_',{'spktm': one_simulation})


# ## One simulation
# one_simulation =  simulation(name_conections='connections_sp_1000_.npz', N=1000, K=25)
# io.savemat('/home/david/Desktop/brian_simulations/single_simulation_1000_',{'spktm': one_simulation})



