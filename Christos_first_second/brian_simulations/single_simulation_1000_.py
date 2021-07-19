#import brian_no_units


from model2 import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3



## One simulation
# one_simulation =  simulation2(loadconnections=False, saveconnections=True, save_name='connections_N1000')
# io.savemat('/home/david/Desktop/brian_simulations/single_simulation',{'spktm': one_simulation})


## One simulation
one_simulation =  simulation(name_conections='connections_N1000.npz', N=1000, K=25)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation',{'spktm': one_simulation})



