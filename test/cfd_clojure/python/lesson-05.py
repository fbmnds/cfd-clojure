from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots

import numpy as np
import matplotlib.pyplot as plt

f = open('test-05.json','w')

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


###variable declarations
nx = 41
ny = 21
nt = 0
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = .2
dt = sigma*dx

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

u = np.ones((ny,nx)) ##create a 1xn vector of 1's
un = np.ones((ny,nx)) ##

f.write("{ \"nx\" : ")
f.write(str(nx))
f.write(", \"dx\" : ")
f.write(str(dx))
f.write(", \"ny\" : ")
f.write(str(ny))
f.write(", \"dy\" : ")
f.write(str(dy))
f.write(", \"dt\" : ")
f.write(str(dt))
f.write(", \"c\" : ")
f.write(str(c))
f.write(", \"sigma\" : ")
f.write(str(sigma))

###Assign initial conditions


u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

print_u (f, ", \"u0\" : [", u, ny, nx, "],\n")


###Plot Initial Condition
fig = plt.figure(figsize=(11,7), dpi=100)          ##the figsize parameter can be used to produce different sized images
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
fig.show()
input(" press ENTER to continue ")



for n in range(nt+1): ##loop across number of time steps
    un[:] = u[:]
    u[1:,1:]=un[1:,1:]-(c*dt/dx*(un[1:,1:]-un[0:-1,1:]))-(c*dt/dy*(un[1:,1:]-un[1:,0:-1]))
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

fig2 = plt.figure(figsize=(11,7), dpi=100)
ax2 = fig2.gca(projection='3d')
surf2 = ax2.plot_surface(X,Y,u[:])
fig2.show()
input(" press ENTER to continue ")

# ... printing u takes time ...
print_u(f, "\n \"u_nt\" : [", u[:],ny, nx, "] } ")

f.close()
