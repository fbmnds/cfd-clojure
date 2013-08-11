from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots

import numpy as np
import matplotlib.pyplot as plt

f = open('test-05.json','w')

###

def print_u (f, label1, u, nx, ny, label2):
    f.write(label1)
    for i in range(ny):
        for j in range(nx):
            f.write(str(u[j,:]))
        f.write("\n")
    f.write(label2)


###variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = .2
dt = sigma*dx

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

u = np.ones((ny,nx)) ##create a 1xn vector of 1's
un = np.ones((ny,nx)) ##


###Assign initial conditions

u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

print_u (f, "{ \"u0\" : [", u, nx, 1, "],\n")


###Plot Initial Condition
fig = plt.figure(figsize=(11,7), dpi=100)          ##the figsize parameter can be used to produce different sized images
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])

u = np.ones((ny,nx))
u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2

for n in range(nt+1): ##loop across number of time steps
    un[:] = u[:]
    u[1:,1:]=un[1:,1:]-(c*dt/dx*(un[1:,1:]-un[0:-1,1:]))-(c*dt/dy*(un[1:,1:]-un[1:,0:-1]))
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

fig = plt.figure(figsize=(11,7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X,Y,u[:])

fig.show()
#input(" press ENTER to continue ")

# ... printing u takes time ...
print_u(f, "\n \"u_nt\" : [", u[:],ny,nx, "] } ")

f.close()
