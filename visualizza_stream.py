import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
nx=64
lx=1
dx=lx/nx
t = FortranFile('temp.bin', 'r')
vey = FortranFile('vey.bin', 'r')
vex = FortranFile('vex.bin', 'r')
datay = vey.read_reals('float')
datax = vex.read_reals('float')
datat = t.read_reals('float')
#vey = datay.reshape(nx, nx)
#vex = datax.reshape(nx, nx)


#vex=np.fromfile("vex.bin","float64",-1)
#vey=np.fromfile("vey.bin","float64",-1)
Y,X=np.mgrid[0:1:dx,0:1:dx]
vex=np.reshape(datax,(nx,nx),order='C')
vey=np.reshape(datay,(nx,nx),order='C')
t=np.reshape(datat,(nx,nx),order='C')
plt.figure(1)
plt.streamplot(X, Y, vey[::-1,::-1],vex[::-1,::-1], density=[4.5, 4.5],linewidth=0.9,arrowsize=0)
plt.figure(2)
plt.contour(X,Y,t,30)
plt.show()
