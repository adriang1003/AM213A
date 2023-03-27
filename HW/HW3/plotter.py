# File: plotter.py
# Author: Adrian Garcia
# Date 02/13/2023
# Purpose: Plot the plane that intersects points A,B,C
import numpy as np 
import matplotlib.pyplot as plt
data = np.genfromtxt('myplot.dat')
x = np.linspace(-10,10,100)
y = np.linspace(-10,10,100)
x,y = np.meshgrid(x,y)
a = data[0]
b = data[1]
c = data[2]

z = -(a/c)*x-(b/c)*y+1/c

fig = plt.figure(figsize = (8,8))
ax = plt.axes(projection = '3d')
surf = ax.plot_surface(x,y,z,alpha=0.6)
ax.scatter(1,2,3,c='r',s=50,label = 'A')
ax.scatter(-3,2,5,c='b',s=50,label = 'B')
ax.scatter(np.pi,np.e,-np.sqrt(2),c='g',s=50,label = 'C')
plt.legend(loc = 'upper right')
ax.view_init(-140,0)
plt.savefig("planeFit.jpg")