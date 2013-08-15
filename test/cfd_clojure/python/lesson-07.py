import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D ##library for 3d projection plots
from matplotlib import cm ##cm = "colormap" for changing the 3d plot color palette

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
nx = 31
ny = 31
nt = 10
nu=.05
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = .25
dt = sigma*dx*dy/nu

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

u = np.ones((ny,nx)) ##create a 1xn vector of 1's
un = np.ones((ny,nx)) ##


f = open(''.join(['test-07-', str(nx), '.json']),'w')

f.write("{ \"nx\" : ")
f.write(str(nx))
f.write(", \"dx\" : ")
f.write(str(dx))
f.write(", \"ny\" : ")
f.write(str(ny))
f.write(", \"dy\" : ")
f.write(str(dy))
f.write(", \"nt\" : ")
f.write(str(nt))
f.write(", \"dt\" : ")
f.write(str(dt))
f.write(", \"nu\" : ")
f.write(str(nu))
f.write(", \"sigma\" : ")
f.write(str(sigma))



###Assign initial conditions

u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

print_u (f, ", \"u0\" : [", u, ny, nx, "],\n")


fig = plt.figure()
ax = fig.gca(projection='3d')
X,Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
plt.show()
ax.set_xlim(1,2)
ax.set_ylim(1,2)
ax.set_zlim(1,2.5)
#ax.zaxis.set_major_locator(LinearLocator(5))


###Run through nt timesteps
def diffuse(nt):
    u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2

    for n in range(nt+1):
        un[:] = u[:]
        u[1:-1,1:-1]=un[1:-1,1:-1]+nu*dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
        u[0,:]=1
        u[-1,:]=1
        u[:,0]=1
        u[:,-1]=1


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)
    ax.set_zlim(1,2.5)
    plt.show()


diffuse(nt)

print_u(f, "\n \"u_nt\" : [", u[:],ny, nx, "] } ")

f.close()
