from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np


###

def print_u (f, label1, u, ny, nx, label2):
    f.write(label1)
    for j in range(ny):
        f.write("[")
        for i in range(nx):
            f.write(str(u[j,i]))
            if i < nx-1:
                f.write(", ")
        f.write("]\n")
    f.write(label2)

def print_v (f, label1, v, n, label2):
    f.write(label1)
    for i in range(n):
        f.write(str(v[i]))
        if i < n-1:
            f.write(", ")
    f.write(label2)

###


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



# Parameters
nx = 50    ## rows
ny = 50    ## cols
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


f = open(''.join(['test-10-', str(nx), '.json']),'w')

f.write("{ \"nx\" : ")
f.write(str(nx))
f.write(", \"dx\" : ")
f.write(str(dx))

f.write(", \"xmax\" : ")
f.write(str(xmax))
f.write(", \"xmin\" : ")
f.write(str(xmin))

f.write(", \"ny\" : ")
f.write(str(ny))
f.write(", \"dy\" : ")
f.write(str(dy))

f.write(", \"ymax\" : ")
f.write(str(ymax))
f.write(", \"ymin\" : ")
f.write(str(ymin))

f.write(", \"nt\" : ")
f.write(str(nt))

print_u (f, ", \"p0\" : [", p, nx, ny, "]\n")
print_u (f, ", \"b\" : [", b, nx, ny, "]\n")


for it in range(nt):

    pd[:][:]=p[:][:]

    p[1:nx-1,1:ny-1] = ( dy**2/(2*(dx**2+dy**2))*(pd[2:nx,1:ny-1]+pd[0:nx-2,1:ny-1]) +
                         dx**2/(2*(dx**2+dy**2))*(pd[1:nx-1,2:ny]+pd[1:nx-1,0:ny-2]) -
                        b[1:nx-1,1:ny-1]*dx**2*dy**2/(2*(dx**2+dy**2)) )

    p[0,:] = p[nx-1,:] = p[:,0] = p[:,ny-1] = 0.0


print_u (f, ", \"p\" : [", p, nx, ny, "] } \n")
f.close()

plot2D(x, y, p)
#input(" press ENTER to continue ")
