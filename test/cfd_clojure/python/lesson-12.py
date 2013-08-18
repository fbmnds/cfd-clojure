from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np


def buildUpB(rho, dt, dx, dy, u, v):
    b = np.zeros_like(u)

    b[1:-1,1:-1]=rho*(1/dt*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx)+(v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))-\
                ((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx))**2-\
                2*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dy)*(v[2:,1:-1]-v[0:-2,1:-1])/(2*dx))-\
                ((v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))**2)

        ####Periodic BC Pressure @ x = 2
    b[-1,1:-1]=rho*(1/dt*((u[0,1:-1]-u[-2,1:-1])/(2*dx)+(v[-1,2:]-v[-1,0:-2])/(2*dy))-\
                ((u[0,1:-1]-u[-2,1:-1])/(2*dx))**2-\
                2*((u[-1,2:]-u[-1,0:-2])/(2*dy)*(v[0,1:-1]-v[-2,1:-1])/(2*dx))-\
                ((v[-1,2:]-v[-1,0:-2])/(2*dy))**2)

        ####Periodic BC Pressure @ x = 0
    b[0,1:-1]=rho*(1/dt*((u[1,1:-1]-u[-1,1:-1])/(2*dx)+(v[0,2:]-v[0,0:-2])/(2*dy))-\
                ((u[1,1:-1]-u[-1,1:-1])/(2*dx))**2-\
                2*((u[0,2:]-u[0,0:-2])/(2*dy)*(v[1,1:-1]-v[-1,1:-1])/(2*dx))-\
                ((v[0,2:]-v[0,0:-2])/(2*dy))**2)

    return b


def presPoissPeriodic(p, dx, dy):
    pn = np.empty_like(p)

    for q in range(nit):
                pn[:]=p[:]
                p[1:-1,1:-1] = ((pn[2:,1:-1]+pn[0:-2,1:-1])*dy**2+(pn[1:-1,2:]+pn[1:-1,0:-2])*dx**2)/\
                        (2*(dx**2+dy**2)) -\
                        dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]

                ####Periodic BC Pressure @ x = 2
                p[-1,1:-1] = ((pn[0,1:-1]+pn[-2,1:-1])*dy**2+(pn[-1,2:]+pn[-1,0:-2])*dx**2)/\
                        (2*(dx**2+dy**2)) -\
                        dx**2*dy**2/(2*(dx**2+dy**2))*b[-1,1:-1]

                ####Periodic BC Pressure @ x = 0
                p[0,1:-1] = ((pn[1,1:-1]+pn[-1,1:-1])*dy**2+(pn[0,2:]+pn[0,0:-2])*dx**2)/\
                        (2*(dx**2+dy**2)) -\
                        dx**2*dy**2/(2*(dx**2+dy**2))*b[0,1:-1]

                ####Wall boundary conditions, pressure
                p[-1,:] =p[-2,:]	##dp/dy = 0 at y = 2
                p[0,:] = p[1,:]		##dp/dy = 0 at y = 0

    return p



##variable declarations
nx = 41
ny = 41
nt = 10
nit=50
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
Y,X = np.meshgrid(y,x)


##physical variables
rho = 1
nu = .1
F = 1
dt = .01

#initial conditions
u = np.zeros((ny,nx)) ##create a XxY vector of 0's
un = np.zeros((ny,nx)) ##create a XxY vector of 0's

v = np.zeros((ny,nx)) ##create a XxY vector of 0's
vn = np.zeros((ny,nx)) ##create a XxY vector of 0's

p = np.ones((ny,nx)) ##create a XxY vector of 0's
pn = np.ones((ny,nx)) ##create a XxY vector of 0's

b = np.zeros((ny,nx))


udiff = 1
stepcount = 0

while udiff > .001:
    un[:] = u[:]
    vn[:] = v[:]

    b = buildUpB(rho, dt, dx, dy, u, v)
    p = presPoissPeriodic(p, dx, dy)

    u[1:-1,1:-1] = un[1:-1,1:-1]-\
                un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
                vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
                dt/(2*rho*dx)*(p[2:,1:-1]-p[0:-2,1:-1])+\
                nu*(dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+\
                dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2]))+F*dt

    v[1:-1,1:-1] = vn[1:-1,1:-1]-\
                un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
                vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
                dt/(2*rho*dy)*(p[1:-1,2:]-p[1:-1,0:-2])+\
                nu*(dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])+\
                (dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])))

        ####Periodic BC u @ x = 2
    u[-1,1:-1] = un[-1,1:-1]-\
                un[-1,1:-1]*dt/dx*(un[-1,1:-1]-un[-2,1:-1])-\
                vn[-1,1:-1]*dt/dy*(un[-1,1:-1]-un[-1,0:-2])-\
                dt/(2*rho*dx)*(p[0,1:-1]-p[-2,1:-1])+\
                nu*(dt/dx**2*(un[0,1:-1]-2*un[-1,1:-1]+un[-2,1:-1])+\
                dt/dy**2*(un[-1,2:]-2*un[-1,1:-1]+un[-1,0:-2]))+F*dt

        ####Periodic BC u @ x = 0
    u[0,1:-1] = un[0,1:-1]-\
                un[0,1:-1]*dt/dx*(un[0,1:-1]-un[-1,1:-1])-\
                vn[0,1:-1]*dt/dy*(un[0,1:-1]-un[0,0:-2])-\
                dt/(2*rho*dx)*(p[1,1:-1]-p[-1,1:-1])+\
                nu*(dt/dx**2*(un[1,1:-1]-2*un[0,1:-1]+un[-1,1:-1])+\
                dt/dy**2*(un[0,2:]-2*un[0,1:-1]+un[0,0:-2]))+F*dt

        ####Periodic BC v @ x = 2
    v[-1,1:-1] = vn[-1,1:-1]-\
                un[-1,1:-1]*dt/dx*(vn[-1,1:-1]-vn[-2,1:-1])-\
                vn[-1,1:-1]*dt/dy*(vn[-1,1:-1]-vn[-1,0:-2])-\
                dt/(2*rho*dy)*(p[-1,2:]-p[-1,0:-2])+\
                nu*(dt/dx**2*(vn[0,1:-1]-2*vn[-1,1:-1]+vn[-2,1:-1])+\
                (dt/dy**2*(vn[-1,2:]-2*vn[-1,1:-1]+vn[-1,0:-2])))

        ####Periodic BC v @ x = 0
    v[0,1:-1] = vn[0,1:-1]-\
                un[0,1:-1]*dt/dx*(vn[0,1:-1]-vn[-1,1:-1])-\
                vn[0,1:-1]*dt/dy*(vn[0,1:-1]-vn[0,0:-2])-\
                dt/(2*rho*dy)*(p[0,2:]-p[0,0:-2])+\
                nu*(dt/dx**2*(vn[1,1:-1]-2*vn[0,1:-1]+vn[-1,1:-1])+\
                (dt/dy**2*(vn[0,2:]-2*vn[0,1:-1]+vn[0,0:-2])))

        ####Wall BC: u,v = 0 @ y = 0,2
    u[:,0] = 0
    u[:,-1] = 0
    v[:,0] = 0
    v[:,-1]=0

    udiff = (np.sum(u)-np.sum(un))/np.sum(u)
    stepcount += 1


print( stepcount )

fig = plt.figure(figsize = (11,7), dpi=100)
plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
fig.show()
input(" press ENTER to continue ")


fig = plt.figure(figsize = (11,7), dpi=100)
plt.quiver(X, Y, u, v)
fig.show()
input(" press ENTER to continue ")
