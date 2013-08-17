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
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.view_init(30,225)
    fig.show()



def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    pn = np.empty_like(p)

    while l1norm > l1norm_target:
        pn[:] = p[:]
        p[1:-1,1:-1] = (dy**2*(pn[2:,1:-1]+pn[0:-2,1:-1])+dx**2*(pn[1:-1,2:]+pn[1:-1,0:-2]))/(2*(dx**2+dy**2))
        p[0,0] = (dy**2*(pn[1,0]+pn[-1,0])+dx**2*(pn[0,1]+pn[0,-1]))/(2*(dx**2+dy**2))
        p[-1,-1] = (dy**2*(pn[0,-1]+pn[-2,-1])+dx**2*(pn[-1,0]+pn[-1,-2]))/(2*(dx**2+dy**2))

        p[:,0] = 0		##p = 0 @ x = 0
        p[:,-1] = y		##p = y @ x = 2
        p[0,:] = p[1,:]		##dp/dy = 0 @ y = 0
        p[-1,:] = p[-2,:]	##dp/dy = 0 @ y = 1
        l1norm = (np.sum(np.abs(p[:])-np.abs(pn[:])))/np.sum(np.abs(pn[:]))

    return p




##variable declarations
nx = 31
ny = 31
c = 1            ## unused
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
eps = .01

f = open(''.join(['test-09-', str(nx), '.json']),'w')

f.write("{ \"nx\" : ")
f.write(str(nx))
f.write(", \"dx\" : ")
f.write(str(dx))
f.write(", \"ny\" : ")
f.write(str(ny))
f.write(", \"dy\" : ")
f.write(str(dy))
f.write(", \"eps\" : ")
f.write(str(eps))

f.write(", \"l1norm\" : 1.")  ## hidden parameter



##initial conditions
p = np.zeros((ny,nx)) ##create a XxY vector of 0's

##plotting aids
x = np.linspace(0,2,nx)
y = np.linspace(0,1,ny)

##boundary conditions
p[:,0] = 0		##p = 0 @ x = 0
p[:,-1] = y		##p = y @ x = 2
p[0,:] = p[1,:]		##dp/dy = 0 @ y = 0
p[-1,:] = p[-2,:]	##dp/dy = 0 @ y = 1


plot2D(x, y, p)

print_u (f, ", \"p0\" : [", p, ny, nx, "]\n")
print_v (f, ", \"y0\" : [", y, ny, "]\n")

p = laplace2d(p, y, dx, dy, eps)

print_u (f, ", \"p\" : [", p, ny, nx, "] \n")
print_v (f, ", \"y\" : [", y, ny, "] } \n")    ## redundant check
f.close()

plot2D(x, y, p)

#input(" press ENTER to continue ")
