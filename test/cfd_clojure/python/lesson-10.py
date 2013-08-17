from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

# Parameters
nx = 50
ny = 50
nt  = 100
xmin = 0.
xmax = 2.
ymin = 0.
ymax = 1.

dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

# Initialization
p  = np.zeros((nx,ny))
pd = np.zeros((nx,ny))
b  = np.zeros((nx,ny))
x  = np.linspace(xmin,xmax,nx)
y  = np.linspace(xmin,xmax,ny)

# Source
b[nx/4][ny/4]  = 100
b[3*nx/4][3*ny/4] = -100


for it in range(nt):

    pd[:][:]=p[:][:]

    p[1:nx-1,1:ny-1] = ( dy**2/(2*(dx**2+dy**2))*(pd[2:nx,1:ny-1]+pd[0:nx-2,1:ny-1]) +
                         dx**2/(2*(dx**2+dy**2))*(pd[1:nx-1,2:ny]+pd[1:nx-1,0:ny-2]) -
                        b[1:nx-1,1:ny-1]*dx**2*dy**2/(2*(dx**2+dy**2)) )

    p[0,:] = p[nx-1,:] = p[:,0] = p[:,ny-1] = 0.0


def plot2D(x, y, p):
    fig = plt.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface( X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False )
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.view_init(30,225)
    fig.show()


plot2D(x, y, p)
input(" press ENTER to continue ")
