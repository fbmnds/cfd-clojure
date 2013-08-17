from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np



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
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)


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


p = laplace2d(p, y, dx, dy, .01)


plot2D(x, y, p)

input(" press ENTER to continue ")
