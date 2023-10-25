import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
nx=128
lx=1
dx=lx/nx

datax=np.fromfile("vex.bin","float64",-1,offset=0)
datay=np.fromfile("vey.bin","float64",-1,offset=0)
datat=np.fromfile("temp.bin","float64",-1,offset=0)

Y,X=np.mgrid[0:1:dx,0:1:dx]
vex=np.reshape(datax,(nx,nx),order='C')
vey=np.reshape(datay,(nx,nx),order='C')
t=np.reshape(datat,(nx+2,nx+2),order='C')
plt.figure(1)
plt.streamplot(X, Y, vex[::1,::1],vey[::1,::1], density=[4.5, 4.5],linewidth=0.9,arrowsize=0.5)
plt.figure(2)
plt.contour(X,Y,t[1:-1,1:-1],30)
plt.show()
