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
one_simulation =  simulation(name_conections='connections_sp_30000.npz', N=30000, extE=2)
io.savemat('/home/david/Desktop/brian_simulations/single_simulation_ext2',{'spktm': one_simulation})




### Multiple simulations in paralel
extEs = list(np.linspace(0,2,10))
results = Parallel(n_jobs = numcores)(delayed(simulation)(extE=extE, name_conections='connections_sp_30000.npz', N=30000)  for extE in extEs)    

io.savemat('/home/david/Desktop/brian_simulations/results_simulations',{'extEs':extEs, 'spktm': results})

