# -*- coding: utf-8 -*-
"""
@author: David Bestue
"""


import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
import pandas as pd
import os
import statsmodels.api as sm
from joblib import Parallel, delayed
import multiprocessing
from scipy.optimize import curve_fit
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import KFold
import random


numcores = multiprocessing.cpu_count() - 2
numcores


path_CTD = '/home/david/Desktop/IDIBAPS/Gottlib'


##functions
def circdist(a1,a2):
    ## Returns the minimal distance in angles between to angles 
    op1=abs(a2-a1)
    angs=[a1,a2]
    op2=min(angs)+(360-max(angs))
    options=[op1,op2]
    return min(options)


####
def decoder_cv(df_train, df_test, splits=100, percentage_training=0.8):
    #### Input : dataframe with three columns: (spikes, behaviour and neuron label)
    df_train.columns=['firing', 'beh', 'trial']
    df_test.columns=['firing', 'beh', 'trial']    
    ######
    ######## Cross validation ########
    ######## different lengths of taringn and testing #####
    errors_splits=[]
    options = df_train.trial.unique()
    for s in range(splits):
        # training
        trials_train = random.sample(list(options), int(percentage_training*len(options)))
        X_train = df_train.loc[df_train['trial'].isin(trials_train)].firing.values
        y_train = df_train.loc[df_train['trial'].isin(trials_train)].beh.values
        X_test = df_test.loc[~df_test['trial'].isin(trials_train)].firing.values ## the trials not used for training and the rest
        y_test = df_test.loc[~df_test['trial'].isin(trials_train)].beh.values
        ######## Trainning #########
        ## X matrix (intercept and spikes)
        X = np.column_stack([np.ones(np.shape(X_train)[0]), X_train])
        ## Y (sinus and cos)
        sinus =np.sin([np.radians(np.array(y_train)[i]) for i in range(0, len(y_train))])
        cosinus = np.cos([np.radians(np.array(y_train)[i]) for i in range(0, len(y_train))])
        Y = np.column_stack([cosinus, sinus])
        ### one OLS for sin and cos: output: beta of intercetp and bea of spikes (two B intercepts and 2 B for spikes )
        model = sm.OLS(Y, X)
        ##train the model
        fit=model.fit()
        ######### Testing ###########
        X = np.column_stack([np.ones(np.shape(X_test)[0]),X_test])
        p = fit.predict(X)
        x = p[:,0]
        y = p[:,1]
        #####
        ##### Error --> take the resulting vector in sin/cos space
        ### from sin and cos get the angle (-pi, pi)
        #pred_angle = [ np.degrees(np.arctan2(y[i], x[i]) + np.pi) for i in range(0, len(y))]
        pred_angle = [ np.degrees(np.arctan2(y[i], x[i])) for i in range(0, len(y))]
        for i in range(0, len(pred_angle)):
            if pred_angle[i]<0:
                pred_angle[i]=360+pred_angle[i]
        ##
        #
        #print(beh_test)
        error_trial=[ circdist(y_test[i], pred_angle[i]) for i in range(0, len(pred_angle))]
        mean_error = np.round(np.mean(error_trial),2)
        errors_splits.append(mean_error)
    ##
    ## mean of all the splits
    mean_error = round(np.mean(errors_splits),2)
    return mean_error
    

#####################
#####################
##################### Cross temporal decoding of all the condition (put together correct and incorrect)
#####################
#####################
#####################

files_ = ['pfc_100.xlsx', 'pfc_200.xlsx','pfc_300.xlsx', 'pfc_900.xlsx', 
'lip_100.xlsx', 'lip_200.xlsx','lip_300.xlsx', 'lip_900.xlsx']

save_files = ['hm_pfc_100.xlsx', 'hm_pfc_200.xlsx', 'hm_pfc_300.xlsx', 'hm_pfc_900.xlsx', 
'hm_lip_100.xlsx', 'hm_lip_200.xlsx', 'hm_lip_300.xlsx', 'hm_lip_900.xlsx']


for idx, file in enumerate(files_):
    print(idx,file)
    condition = pd.read_excel( os.path.join(path_CTD, file))
    ### condition = condition.loc[condition['performance_code']==1] ##get just the correct ones
    ############## 
    heatmaps_=[]
    for n_neuron,neuron_ in enumerate(condition.neuron.unique()):
        print(n_neuron, neuron_)
        ## get the neuron
        dfN = condition.loc[condition['neuron']==neuron_]
        ## column of times centered to stim onset
        dfN['times_centered'] = dfN['times'] - dfN['fixationtime']
        #
        #cross-decoding for 1 neuron
        ### empty matrix to append the cross temporal decoding
        all_times = np.arange(-500,2100, 100)
        train_test = np.empty( (len(all_times), len(all_times) ) )
        train_test[:] = np.nan
        #
        ### times I am interested in (no need to get more times)
        dfN = dfN.loc[(dfN['times_centered']>=-500) & (dfN['times_centered']<=2000) ]
        list_times = dfN.times_centered.unique()
        list_times_sorted =np.sort(list_times)
        #list_times_sorted
        #
        for training_time in list_times_sorted:
            #print(training_time)
            dfn_train = dfN.loc[(dfN['times_centered']==training_time), ['firing', 'target_angle', 'trial']]         
            dfn_train = dfn_train.loc[dfn_train.iloc[:,0]<9999] ###Take off nans
            paralel_train = [dfn_train for i in range(len(list_times_sorted))]
            #
            paralel_test = []
            for times_testing in list_times_sorted:
                dfn_test =dfN.loc[(dfN['times_centered']==times_testing), ['firing', 'target_angle', 'trial']]
                dfn_test = dfn_test.loc[dfn_test.iloc[:,0]<9999] ###Take off nans
                paralel_test.append(dfn_test)
            #
            cross_temp = Parallel(n_jobs = numcores)(delayed(decoder_cv)(training, testing)  for training, testing in zip(paralel_train, paralel_test) )   #### reconstruction standard (paralel)
            #
            idx_append_row = np.where(training_time == all_times )[0][0]
            idx_append_col1 = np.where(list_times_sorted[0] == all_times)[0][0]
            idx_append_col2 = np.where(list_times_sorted[-1] == all_times)[0][0]
            train_test[idx_append_row, idx_append_col1:idx_append_col2+1] = cross_temp
        ##
        heatmaps_.append(train_test)
    ###
    ###
    path_save = os.path.join(path_CTD, save_files[idx])
    #
    writer = pd.ExcelWriter(path_save)
    for idx, neuron_name in  enumerate(condition.neuron.unique()):
        hm_ = pd.DataFrame(heatmaps_[idx])
        hm_.to_excel(writer, sheet_name=str(neuron_name)) #each dataframe in a excel sheet
    #
    writer.save()   #save reconstructions (heatmaps)
#


#