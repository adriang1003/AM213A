# File: plotter.py
# Author: Adrian Garcia
# Date 02/27/2023
# Purpose: Plots the fitted curve to the given data points in 'atkinson.dat' for degree 3 and 5 polynomials using Cholesky and QR decomposition
# Run Command: Type 'python3 plotter.py' to command window
## Presets
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
def Degree3(constants,x):
    return constants[0]+constants[1]*x+constants[2]*x**2+constants[3]*x**3
def Degree5(constants,x):
    return constants[0]+constants[1]*x+constants[2]*x**2+constants[3]*x**3+constants[4]*x**4+constants[5]*x**5
plt.rcParams["figure.figsize"] = (10, 8)
plt.rcParams.update({'font.size': 12})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXgeneral'
## Load data
data = np.loadtxt('atkinson.dat')
Cho3 = np.loadtxt('ChoDegree3.dat')
Cho5 = np.loadtxt('ChoDegree5.dat')
QR3 = np.loadtxt('QRDegree3.dat')
QR5 = np.loadtxt('QRDegree5.dat')
## Calculate 2-norm error between the fitted curve and the data (Cholesky Decomposition)
print('(Cholesky: Degree 3) Error: ||x||_2 =',LA.norm(Degree3(Cho3,data[1:,0])-data[1:,1]))
print('(Cholesky: Degree 5) Error: ||x||_2 =',LA.norm(Degree5(Cho5,data[1:,0])-data[1:,1]))
## Calculate 2-norm error between the fitted curve and the data (QR Decomposition)
print('(QR: Degree 3) Error: ||x||_2 =',LA.norm(Degree3(QR3,data[1:,0])-data[1:,1]))
print('(QR: Degree 5) Error: ||x||_2 =',LA.norm(Degree5(QR5,data[1:,0])-data[1:,1]))
## Plot data (Cholesky Decomposition)
fig = plt.figure()
plt.plot(data[1:,0],data[1:,1],'o',label = 'Data Points')
plt.plot(data[1:,0],Degree3(Cho3,data[1:,0]),'red',label = '3rd Degree Polynomial')
plt.plot(data[1:,0],Degree5(Cho5,data[1:,0]),'blue',label = '5th Degree Polynomial')
## Plot config (Cholesky Decomposition)
plt.title('Data Fit (Cholesky)', fontsize = 18)
plt.xlabel('$x$', fontsize = 14)
plt.ylabel('$y$', fontsize = 14)
plt.legend(loc = 'upper left')
plt.grid(True)
plt.savefig('ChoPolyFit.jpg')
## Plot data (QR Decomposition)
fig = plt.figure()
plt.plot(data[1:,0],data[1:,1],'o',label = 'Data Points')
plt.plot(data[1:,0],Degree3(QR3,data[1:,0]),'red',label = '3rd Degree Polynomial')
plt.plot(data[1:,0],Degree5(QR5,data[1:,0]),'blue',label = '5th Degree Polynomial')
## Plot config (QR Decomposition)
plt.title('Data Fit (QR)', fontsize = 18)
plt.xlabel('$x$', fontsize = 14)
plt.ylabel('$y$', fontsize = 14)
plt.legend(loc = 'upper left')
plt.grid(True)
plt.savefig('QRPolyFit.jpg')
plt.show()