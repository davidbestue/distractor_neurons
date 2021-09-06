# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:05:33 2019

@author: David
"""

import scipy.io as io
import pickle
from scipy.ndimage import gaussian_filter
from scipy import misc
import os
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from joblib import Parallel, delayed
import multiprocessing


numcores = multiprocessing.cpu_count() 
if numcores>20:
    numcores=numcores-5
if numcores<10:
    numcores=numcores-3



def circ_dist(a1,a2):
    ## Returns the minimal distance in angles between to angles 
    op1=abs(a2-a1)
    angs=[a1,a2]
    op2=min(angs)+(360-max(angs))
    options=[op1,op2]
    return min(options)


def decode_rE(rE, a_ini=0, a_fin=360, N=512):
    #Population vector for a given rE
    # return ( angle in radians, absolut angle in radians, abs angle in degrees )
    N=len(rE)
    Angles = np.linspace(a_ini, a_fin, N) 
    angles=np.radians(Angles)
    rE = np.reshape(rE, (1,N))
    R = np.sum(np.dot(rE,np.exp(1j*angles)))/np.sum(rE)
    
    angle_decoded = np.degrees(np.angle(R))
    if angle_decoded<0:
        angle_decoded = 360+angle_decoded
    
    return angle_decoded


def err_deg(a1,ref):
    ### Calculate the error ref-a1 in an efficient way in the circular space
    ### "going ccw from the reference (ref) is negative" - "going cw from the reference (ref) is positive"
    ### it uses complex numbers!
    ### Input in degrees (0-360)
    a1=np.radians(a1)
    ref=np.radians(ref)
    err = np.angle(np.exp(1j*ref)/np.exp(1j*(a1) ), deg=True) 
    err=round(err, 2)
    return err


###


path_ =  '/home/david/Desktop/brian_simulations_albert/simulations8'
all_= os.listdir(path_)



N0 = 20000
time_s =7000 #(ms)
N=0.8*N0 #(el 80% son excitadoras)
rounding = 2 ##round the timing
w=1000 #(ms)

pos_stim=[]
Iexts = []
firings_wind = []

for sim_ in np.arange(0,250,1): ##simulations used (all or a fraction e.j. np.arange(0,500,100)) )
    #
    print(sim_)
    simx = io.loadmat(path_ + '/' + all_[sim_])
    ####
    #### save position and I0E of each simulation
    ####
    pos_stim.append(360*simx['pos_stim'][0][0])
    Iexts.append(simx['IEext'][0][0])

    spikes = simx['spktm']
    
    Matrix_spikes=np.zeros([int(0.8*N0), time_s])

    neurons_ = np.array([int(spikes[0][x]) for x in range(len(spikes[0]))])

    times_ = np.array([spikes[1][x] for x in range(len(spikes[1]))])
    times_ = times_ - min(times_)
    times_ = np.array([int(times_[x]*1000) for x in range(len(times_))])

    # dimensions = neurons, time_s
    for t, n in zip(times_, neurons_ ):
        Matrix_spikes[n,t]=1
    ##
    ####
    #### STEP 2: calculate firing of each neuron in windows of 100ms (10)
    ####
    t1s = np.arange(0, time_s ,w)
    t2s = np.arange(w, time_s+w, w)
    ##
    fr_time = []
    for N in range(np.shape(Matrix_spikes)[0]):
        neuron_fr = []
        for i in range(len(t1s)):
            neuron_fr.append(Matrix_spikes[N, t1s[i]:t2s[i]].sum()/ float(1000/w) ) 
        #
        fr_time.append(neuron_fr)
    ###
    fr_time=np.array(fr_time)
    ####
    #### STEP 3: save the new matrix, this time with the dimension (neuron, time) ##time dimension is 700/w = 70
    ####
    firings_wind.append(fr_time)
    
    
##
IExts = np.array(Iexts) 
Positions=np.array(pos_stim)
firings_wind = np.array(firings_wind)



##########################
##########################

IExts = np.array(Iexts) 
Positions=np.array(pos_stim)
firings_wind = np.array(firings_wind)


Neurons_  = np.arange(0,16000,1000) ## np.arange(0,16000,1) 
Windows_ =np.arange(0,7000,w)  



Iext_ = [0, 0.5]


lim_RF = 45/4 ##limit to consider a position is inside the RF of the neuron
base_fpr = np.linspace(0, 1, 101)

number_Iexs = len(Iext_)
number_neurons=len(Neurons_)
number_windows = len(Windows_)


auc_ = np.zeros((number_Iexs, number_windows, number_neurons))
#tprs_ =np.zeros((number_Iexs, number_windows))


for idx_Iext, IEXT in enumerate(Iext_):
    for idx_wind, wind in enumerate(Windows_):
        print(IEXT, wind)
        for idx_neuron, Neuron in enumerate(Neurons_):
            ###
            ### Get the positions of the appropiate Iext
            nx_positions=Positions[IExts==IEXT]
            ###
            ### Get the simulations of the appropiate Iext
            nx_rates_iex = firings_wind[IExts==IEXT]
            ###
            ### Get the firing of the neuron at a certain window in each simulation
            nx_rates = np.array([nx_rates_iex[n][Neuron, idx_wind] for n in range(len(nx_rates_iex))])
            ###
            dfx = pd.DataFrame({'position':nx_positions, 'rate':nx_rates, 'Neuron':Neuron, 'wind':wind})
            ###
            ### RF center of each neuro
            index_max_rate = np.where(dfx['rate']==dfx['rate'].max())[0][0]
            RF_center = dfx['position'].iloc[index_max_rate]
            dfx['RF_center'] = RF_center
            in_out_rf = []
            for p in range(len(dfx)):
                dist_ = circ_dist(RF_center, dfx.position.iloc[p])
                if dist_>lim_RF:
                    in_out_rf.append(5) #outside RF
                else:
                    in_out_rf.append(1) #inside RF
            #####
            dfx['in_out_rf']=in_out_rf
            ###
            ## appropiate shape for the classifier
            y = label_binarize(dfx['in_out_rf'].values, classes=[5,1]) #matrix (1,0,0,...,0)
            y=y.ravel()
            ######
            if len(np.unique(y))==1:
                mean_tprs = np.empty(np.shape(base_fpr))*np.nan
            else:
                fpr, tpr, _ = roc_curve(y, dfx['rate'].values) #HERE: compute ROC with raw firing rate
                tpr = np.interp(base_fpr, fpr, tpr)
                tpr[0] = 0.0
                mean_tprs = tpr
            
            #tprs_[idx_Iext, idx_wind] = mean_tprs
            ##
            auc_[idx_Iext, idx_wind, idx_neuron] = auc(base_fpr, mean_tprs)
            ##




### Mean AUC of all neurons in each Iex and window

Mean_auc = []
exts=[]
times_=[]

for idx_Iext, IEXT in enumerate(Iext_):
    for idx_wind, wind in enumerate(Windows_):
        auc_wind = auc_[idx_Iext, idx_wind, :] 
        ind = np.isnan(auc_wind)
        ##
        Mean_auc.append(np.mean(auc_wind[~ind]) )
        exts.append(IEXT)
        times_.append(wind)
        
##

df_mean_auc = pd.DataFrame({'AUC':Mean_auc, 'Iext':exts, 'time':times_})
df_mean_auc[' '] = df_mean_auc['Iext'].replace([0, 0.5], ['OFF', 'ON'])


df_mean_auc.to_excel('/home/david/Desktop/brian_simulations_albert/auc_time_1000_newregime.xlsx')


