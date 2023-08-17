#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import matplotlib.pyplot as plt
import numba as nb


# In[10]:


def e_integrand(zeta, m_e, T_values):
    integrand = np.zeros((len(T_values), len(zeta)))
    
    for i in range(len(T_values)):
        T = T_values[i]
        for j in range(len(zeta)):
            integrand[i, j] = (zeta[j]**2 * np.sqrt(zeta[j]**2 + (m_e / T)**2)) / (np.exp(np.sqrt(zeta[j]**2 + (m_e / T)**2)) + 1) * np.exp(zeta[j])
    
    return integrand

x_gauss,w_gauss = np.polynomial.laguerre.laggauss(40)


# This function below allows us to calcualte dQ/dt and the enegery density in neutrinos (pv) from the mass file. Takes the name of the mass file as an input and the name of the saved file you would like to be created.

# In[11]:


def imp_values(filename,saved_file):
    actual_data= np.load(filename, allow_pickle=True)
    f_array = actual_data['fe']
    e_array = actual_data['e']
    temp_cm = actual_data['Tcm']
    m_s = actual_data['mass']#in MeV
    t = actual_data['time'] #inverse MeV
    a = actual_data['scalefactors']
    life = actual_data['lifetime']/(6.58*10**-25)*1/1000 #seconds need to be converted to inverse MeV
    temp = actual_data['temp'] #units of MeV
    decay = actual_data['decayrate']
    coll = actual_data['collisionrate']
    p_n = actual_data['p_n_rate']
    n_p = actual_data['n_p_rate']
    mix = actual_data['mixing']
    T_0 = temp[0]
    roww = len(temp)
    
    
    eta_0 = 6.1e-10
    fail_safe = 1
               
               
    D = (10.75/61.75)
    Rz_3 = 1.202056901178332
    n_s = D*(3*Rz_3/(2*np.pi**2))*(temp_cm)**3*np.exp(-t/life)
    ns_initial = n_s[0]
    
    
    m_pl = (1.2*10**19)*1000 #planck mass in MeV
    m_e = 0.510 #electron mass in MeV
    
    
    p_y = (np.pi**2/15)*(temp)**4
    p_vs = m_s*n_s

    results = []
    for e,f in zip(e_array,f_array):
        result = np.trapz(e**3*f,e)
        results.append(result)
        
    p_v = 6*((temp_cm)**4/(2*np.pi**2))*results #in units of MeV^4
    
    
    integrand_vals = np.zeros((len(e_integrand(x_gauss,0.510,temp)),len(w_gauss)))
    vals = e_integrand(x_gauss,m_e,temp)
    
    for i in range(len(integrand_vals)):
        V = vals[i]
        for j in range(len(w_gauss)):
            integrand_vals[i,j] = (2*(temp[i])**4/(np.pi)**2)*V[j]*w_gauss[j]
            
    sums_array = []
    for subarray in integrand_vals:
        subarray_sum = 0
        for i in subarray:
            subarray_sum += i
        sums_array.append(subarray_sum)
   
    p_e = sums_array
    p = p_e+p_y+p_v+p_vs
    
    da = np.sqrt(8*np.pi/(3*m_pl**2))*a/(p)**(-1/2)
    
    results_1 = []
    for e,dec,col in zip(e_array,decay,coll):
        result_1 = np.trapz(e**3*(dec+col),e)
        results_1.append(result_1)
    
    dQ = (m_s*n_s*a**3/life)-(temp_cm**4*a**3/(2*np.pi**2))*da*results_1 #in units of MeV^5
    
    dqdt = dQ/a**3 #needed to match Hannah's values 
    
    
   
    np.savez(saved_file, T = temp, Tcm = temp_cm, dQdt = dqdt, pn = p_n, np = n_p, rhonu = p_v, T_initial = T_0, eta_initial = eta_0, failsafe = fail_safe, ns_0 = ns_initial, vs_mass = m_s, mix_ang = mix, row = roww, pe_dens = p_e, scale = a)

    
    
    






