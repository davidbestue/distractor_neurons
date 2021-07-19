#import brian_no_units


from model import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3



## First simulation
one_simulation =  simulation2(N=1000, loadconnections=False, saveconnections=True)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation',{'spktm': one_simulation})


# ## One simulation
# one_simulation =  simulation2(N=80000)
# io.savemat('/home/david/Desktop/brian_simulations/single_simulation',{'spktm': one_simulation})


