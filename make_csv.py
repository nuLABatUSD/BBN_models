import numpy as np
import os
from scipy.interpolate import CubicSpline
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def delete_file(file):
    check_file = Path(file)

    if check_file.is_file():
        os.remove(check_file) 

def delete_csv_files(folder_name):
    delete_file(folder_name + "T_dqdt.csv")
    delete_file(folder_name + "a_dqdt.csv")
    delete_file(folder_name + "b_dqdt.csv")
    delete_file(folder_name + "c_dqdt.csv")
    delete_file(folder_name + "d_dqdt.csv")

    delete_file(folder_name + "T_pn.csv")
    delete_file(folder_name + "a_pn.csv")
    delete_file(folder_name + "b_pn.csv")
    delete_file(folder_name + "c_pn.csv")
    delete_file(folder_name + "d_pn.csv")

    delete_file(folder_name + "T_np.csv")
    delete_file(folder_name + "a_np.csv")
    delete_file(folder_name + "b_np.csv")
    delete_file(folder_name + "c_np.csv")
    delete_file(folder_name + "d_np.csv")

    delete_file(folder_name + "T_rhonu.csv")
    delete_file(folder_name + "a_rhonu.csv")
    delete_file(folder_name + "b_rhonu.csv")
    delete_file(folder_name + "c_rhonu.csv")
    delete_file(folder_name + "d_rhonu.csv")
    
def generate_csv(npz_file, folder_name = "CSVFiles", save = True):
    load_data = np.load(npz_file)
    
    if folder_name[-1] != '/':
        folder_name = folder_name + '/'

    T = load_data['T']
    Tcm = load_data['Tcm']
    dqdt = load_data['dQdt']
    
    pn = load_data['pn']
    NP = load_data['np']
    rhonu = load_data['rhonu']

    T_reversed = T[::-1]
    Tcm_reversed = Tcm[::-1]
    dqdt_reversed = dqdt[::-1]
    pn_reversed = pn[::-1]
    np_reversed = NP[::-1]
    rhonu_reversed = rhonu[::-1]

    T_fit_pn = CubicSpline(np.log(T_reversed[1:]), np.log(pn_reversed[1:])) #cubic spline fit 
    T_fit_np = CubicSpline(np.log(T_reversed), np.log(np_reversed)) #cubic spline fit 
    
    T_fit_dqdt = CubicSpline(T_reversed[:-1], dqdt_reversed[:-1]) #cubic spline fit 
    Tcm_fit_rhonu = CubicSpline(Tcm_reversed, rhonu_reversed) #cubic spline fit 
    
    #determine coefficients
    d_dqdt = T_fit_dqdt(T)
    c_dqdt = T_fit_dqdt(T,1)
    b_dqdt = T_fit_dqdt(T,2)*.5
    a_dqdt = T_fit_dqdt(T,3)*(1/6)
            
    d_pn = T_fit_pn(np.log(T))
    c_pn = T_fit_pn(np.log(T),1)
    b_pn = T_fit_pn(np.log(T),2)*.5
    a_pn = T_fit_pn(np.log(T),3)*(1/6)

    d_np = T_fit_np(np.log(T))
    c_np = T_fit_np(np.log(T),1)
    b_np = T_fit_np(np.log(T),2)*.5
    a_np = T_fit_np(np.log(T),3)*(1/6)

    d_rhonu = Tcm_fit_rhonu(Tcm)
    c_rhonu = Tcm_fit_rhonu(Tcm,1)
    b_rhonu = Tcm_fit_rhonu(Tcm,2)*.5
    a_rhonu = Tcm_fit_rhonu(Tcm,3)*(1/6)
    
    try:
        os.mkdir(folder_name)
    except OSError as e:
        print("Directory exists")
        
    delete_csv_files(folder_name)
    
    T_dqdt_csv = pd.DataFrame(T).to_csv(folder_name + 'T_dqdt.csv')
    a_dqdt_csv = pd.DataFrame(a_dqdt).to_csv(folder_name + 'a_dqdt.csv')
    b_dqdt_csv = pd.DataFrame(b_dqdt).to_csv(folder_name + 'b_dqdt.csv')
    c_dqdt_csv = pd.DataFrame(c_dqdt).to_csv(folder_name + 'c_dqdt.csv')
    d_dqdt_csv = pd.DataFrame(d_dqdt).to_csv(folder_name + 'd_dqdt.csv')

    T_pn_csv = pd.DataFrame(T).to_csv(folder_name + 'T_pn.csv')
    a_pn_csv = pd.DataFrame(a_pn).to_csv(folder_name + 'a_pn.csv')
    b_pn_csv = pd.DataFrame(b_pn).to_csv(folder_name + 'b_pn.csv')
    c_pn_csv = pd.DataFrame(c_pn).to_csv(folder_name + 'c_pn.csv')
    d_pn_csv = pd.DataFrame(d_pn).to_csv(folder_name + 'd_pn.csv')

    T_np_csv = pd.DataFrame(T).to_csv(folder_name + 'T_np.csv')
    a_np_csv = pd.DataFrame(a_np).to_csv(folder_name + 'a_np.csv')
    b_np_csv = pd.DataFrame(b_np).to_csv(folder_name + 'b_np.csv')
    c_np_csv = pd.DataFrame(c_np).to_csv(folder_name + 'c_np.csv')
    d_np_csv = pd.DataFrame(d_np).to_csv(folder_name + 'd_np.csv')

    Tcm_rhonu_csv = pd.DataFrame(Tcm).to_csv(folder_name + 'T_rhonu.csv')
    a_rhonu_csv = pd.DataFrame(a_rhonu).to_csv(folder_name + 'a_rhonu.csv')
    b_rhonu_csv = pd.DataFrame(b_rhonu).to_csv(folder_name + 'b_rhonu.csv')
    c_rhonu_csv = pd.DataFrame(c_rhonu).to_csv(folder_name + 'c_rhonu.csv')
    d_rhonu_csv = pd.DataFrame(d_rhonu).to_csv(folder_name + 'd_rhonu.csv')
    
    if save == True:
        np.savez(folder_name + "csvFiles.npz", *[T, Tcm, a_dqdt, b_dqdt, c_dqdt, d_dqdt, a_pn, b_pn, c_pn, d_pn, a_np, b_np, c_np, d_np, a_rhonu, b_rhonu, c_rhonu, d_rhonu])
    
    print(".csv files generated in " + folder_name)
    return 