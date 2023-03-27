# File: plotter.py
# Author: Adrian Garcia
# Date 03/18/2023
## Presets
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXgeneral'
index = [10,100,1000]
filename1 = ['GJError10.dat','GJError100.dat','GJError1000.dat']
filename2 = ['GSError10.dat','GSError100.dat','GSError1000.dat']
for i in range(len(filename1)):
    ## Load data
    data1 = np.loadtxt(filename1[i],dtype = np.float32)
    data2 = np.loadtxt(filename2[i],dtype = np.float32)
    ## Figure config
    fig = plt.figure(i+1)
    ## Plot data
    plt.plot(data1[:,1],data1[:,2],'o',label = 'Gauss-Jacobi')
    plt.plot(data2[:,1],data2[:,2],'o',label = 'Gauss-Seidel')
    ## Plot config
    plt.yscale('log')
    plt.legend()
    plt.suptitle('$D$ = %d'%index[i],fontsize = 18)
    plt.savefig('Convergence%d.jpg'%index[i],bbox_inches='tight',dpi=300)