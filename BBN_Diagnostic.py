#!/usr/bin/env python
# coding: utf-8

# In[22]:


import numpy as np
import matplotlib.pyplot as plt
import numba as nb
import subprocess
import ipywidgets as widgets
import pandas as pd
import make_csv
import BBN_Values as BBN_V


# The function below takes a mass file as an input as well as a file-name that you would like to be saved. It will be saved as a .npz file and should be the input to generate_2 and generate_3

# In[23]:


def generate_1(filename,saved_file): 
    BBN_V.imp_values(filename,saved_file)


# This function takes saved_file and generates all the CSV files. In order to make sure the function has ran correctly, go to downloads-->BBN_model-->alterbbn_v2.2-->CSVFiles. All the files here should have been updated seconds ago if the function is running correctly. 

# In[24]:


def generate_2(saved_file):
    make_csv.generate_csv(saved_file, folder_name="alterbbn_v2.2/CSVFiles/")


# This function takes saved_file as in input as well. It then runs the file alter_vs.c to generate a new evolution_vs.out file. In order to make sure this function has run correctly go to downloads-->BBN_models. The file called evolution_vs.out should have been updated seconds ago if the function was running correctly with no errors. 

# In[25]:


def generate_3(saved_file):
    actual_data= np.load(saved_file, allow_pickle=True)
    T_init = actual_data['T_initial']
    eta_init = actual_data['eta_initial']
    fail_safe = actual_data['failsafe']
    ns_init = actual_data['ns_0']
    m_s = actual_data['vs_mass']
    mix = actual_data['mix_ang']
    row = actual_data['row']

    subprocess.run(['./alterbbn_v2.2/alter_vs.x', str(T_init),str(eta_init),str(fail_safe),str(row),str(ns_init),str(m_s),str(mix)], capture_output = True, text = False)


# generates diagnostics as a test 

# In[26]:


def generate_4(filename,saved_file,file_location):
    actual_data= np.load(filename, allow_pickle=True)
    f_array = actual_data['fe']
    e_array = actual_data['e']
    temp_cm = actual_data['Tcm']
    m_s = actual_data['mass']#in MeV
    t = actual_data['time']*1000*6.582*10**-25#conversion to seconds for plots 
    a = actual_data['scalefactors']
    life = actual_data['lifetime']/(6.58*10**-25)*1/1000 #seconds need to be converted to inverse MeV
    temp = actual_data['temp'] #units of MeV
    decay = actual_data['decayrate']
    coll = actual_data['collisionrate']
    p_n = actual_data['p_n_rate']
    n_p = actual_data['n_p_rate']
    T_reversed = temp[::-1]

    
    more_data= np.load(saved_file, allow_pickle=True)
    p_e = more_data['pe_dens']
    p_e_reversed = p_e[::-1]

    
    
    data_file = file_location
    time,scale,temperature = np.loadtxt(data_file,delimiter=',',skiprows=1,usecols=(0,1,2),unpack=True)
    
    plt.figure()
    plt.loglog(time,(scale*temperature)/(scale[0]*temperature[0]),color='red')
    plt.loglog(t,(a*temp)/(a[0]*temp[0]),color='red',linestyle='--')
    plt.show()
    
    plt.figure()
    plt.loglog(T_reversed,p_e_reversed/(T_reversed**4),color='red')
    plt.show()
    






