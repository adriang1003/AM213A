# File: plotter.py
# Author: Adrian Garcia
# Date 03/18/2023
## Presets
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXgeneral'
## Plotting SVD Decomp
k = [10,20,40,80,160,320,640,1279]
filename = ['Image_appn_100010.dat','Image_appn_100020.dat','Image_appn_100040.dat','Image_appn_100080.dat',\
            'Image_appn_100160.dat','Image_appn_100320.dat','Image_appn_100640.dat','Image_appn_101279.dat','dog_bw_data.dat']
fig = plt.figure(1)
fig.tight_layout()
for i in range(len(filename)):
    ## Load data
    if filename[i] == 'dog_bw_data.dat':
        data = np.loadtxt(filename[i])
    else:
        data = np.transpose(np.loadtxt(filename[i]))
    ## Plot data
    plt.subplot(3,3,i+1)
    plt.imshow(data)
    ## Plot config
    plt.set_cmap('gist_gray')
    plt.xticks([])
    plt.yticks([])
    if filename[i] == 'dog_bw_data.dat':
        plt.title('Original',fontsize = 10)
    else:
        plt.title('$\sigma_k$ = %d'%k[i],fontsize = 10)
## Plot config (cont.)
plt.suptitle('SVD Decomposition',fontsize = 18)
plt.savefig('SVDDecomp.jpg',bbox_inches='tight',dpi=300)
## Plotting errors
fig = plt.figure(2)
fig.tight_layout()
## Load data
data = np.loadtxt('error.dat')
## Plot data
plt.plot(data[:,0],data[:,1])
## Plot config
plt.grid('True')
plt.title('$\sigma_k$ vs. $E_k$',fontsize = 18)
plt.xlabel('$\sigma_k$',fontsize = 12)
plt.ylabel('$E_k$',fontsize = 12)
plt.savefig('SVDError.jpg',bbox_inches='tight',dpi=300)