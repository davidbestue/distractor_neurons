#import brian_no_units


from model import * 
from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-10
if numcores<10:
    numcores=numcores-3



### One simulation
one_simulation =  simulation(name_conections='connections_sp.npz', N=1000)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation',{'spktm': one_simulation})



### Multiple simulations in paralel
extEs = [0,1,2]
results = Parallel(n_jobs = numcores)(delayed(simulation)(extE=extE, name_conections='connections_sp.npz', N=1000)  for extE in extEs)    

io.savemat('/home/david/Desktop/brian_simulations/results_simulations',{'extEs':extEs, 'spktm': results})


