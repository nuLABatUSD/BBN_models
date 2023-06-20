import subprocess
import numpy as np
import matplotlib.pyplot as plt


def clean_directory(directory, file, if_PC = True):
    # enter your file directory for alterbbn_v2.2 
    #run_directory = 'C:/Users/kathr/Chad_Research_Neutrinos/alterbbn_v2.2/alter_eta_edit.x'
    #run_directory_new = 'C:/Users/kathr/Chad_Research_Neutrinos/alterbbn_v2.2/alter_vs.x'

    #directory = 'C:/Users/kathr/Chad_Research_Neutrinos/alterbbn_v2.2'

    if if_PC or directory != '':
        print(['make', '-C', directory, file])
        ch = subprocess.run(['make','-C', directory, 'clean'], capture_output = True, text = True, shell = if_PC)
        sh = subprocess.run(['make', '-C', directory, file], capture_output = True, text = True, shell = if_PC)
    else:
        ch = subprocess.run(['make', 'clean'], capture_output = True, text = True, shell = False)
        sh = subprocess.run(['make', file], capture_output = True, text = True, shell = False)
    print(ch.stderr)
    print(sh.stderr)

def output_plots_alter_vs(nu_start, nu_end, nu_steps, run_directory_new, init_T, row, ns0, mass, mixangle, if_PC = True):   
    #nu = np.linspace(1e-10, 10e-10, 40)
    nu = np.linspace(nu_start, nu_end, nu_steps)

    Yp_vs = np.zeros(len(nu))
    H2_H_vs = np.zeros(len(nu))
    He3_H_vs = np.zeros(len(nu))
    Li7_H_vs = np.zeros(len(nu))
    Li6_H_vs = np.zeros(len(nu))
    Be7_H_vs = np.zeros(len(nu))

    a = 0
#'14.87365053664758'
    for i in nu:
        d = subprocess.run([run_directory_new, init_T, str(i), '1', row, ns0, mass, mixangle], capture_output = True, text = True, shell = if_PC)
        #print(d)
        #print(d.stderr)
        x = d.stdout
        #print(x)
        s = x.split('\t')
        #print(s)

        Yp_vs[a] = float(s[12])
        H2_H_vs[a] = float(s[13])
        He3_H_vs[a] = float(s[14])
        Li7_H_vs[a] = float(s[15])
        Li6_H_vs[a] = float(s[16])
        #Be7_H[a] = float(s[17])
        a += 1
        
    #print(np.abs(Li7_H_vs))
    
    plt.figure()
    plt.loglog(nu, Li7_H_vs)
    plt.ylabel('Li7/H')
    plt.xlabel('$\eta_0$')
    #plt.ylim(10e-11, 10e-9)
    #plt.xlim(1e-10, 10e-10)
    #plt.savefig('Li7.jpg')
    plt.show()
    
    #print(np.abs(He3_H_vs))
    
    #plt.figure()
    plt.loglog(nu, He3_H_vs)
    plt.ylabel('He3/H')
    plt.xlabel('$\eta_0$')
    #plt.ylim(10e-7, 10e-4)
    #plt.xlim(1e-10, 10e-10)
    #plt.savefig('He3.jpg')
    plt.show()
    
    np.savez("vs.npz", eta = nu, Li = Li7_H_vs, D = H2_H_vs)
    


def output_plots_alter_eta(nu_start, nu_end, nu_steps, run_directory, if_PC = True):
    nu = np.linspace(nu_start, nu_end, nu_steps)
    
    Yp = np.zeros(len(nu))
    H2_H = np.zeros(len(nu))
    He3_H = np.zeros(len(nu))
    Li7_H = np.zeros(len(nu))
    Li6_H = np.zeros(len(nu))
    Be7_H = np.zeros(len(nu))

    a = 0

    for i in nu:
        d = subprocess.run([run_directory, str(i)], capture_output = True, text = True, shell = if_PC)
        #print(d.stderr)
        x = d.stdout
        #print(x)
        s = x.split('\t')
        #print(s)

        Yp[a] = float(s[12])
        H2_H[a] = float(s[13])
        He3_H[a] = float(s[14])
        Li7_H[a] = float(s[15])
        Li6_H[a] = float(s[16])
        #Be7_H[a] = float(s[17])
        a += 1

    plt.figure()
    plt.loglog(nu, Li7_H)
    plt.ylabel('Li7/H')
    plt.xlabel('$\eta_0$')
    plt.ylim(10e-11, 10e-9)
    plt.xlim(1e-10, 10e-10)
    plt.show()

    plt.figure()
    plt.loglog(nu, H2_H)
    plt.ylabel('D/H')
    plt.xlabel('$\eta_0$')
    plt.ylim(1e-5, 1e-3)
    plt.xlim(1e-10, 10e-10)
    plt.show()
        
    np.savez("eta.npz", eta = nu, Li = Li7_H, D = H2_H)
