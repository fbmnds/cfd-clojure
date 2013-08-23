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

def print_params(f,nx,dx,ny,dy,nt,dt,rho,nu,nit):
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

    f.write(", \"rho\" : ")
    f.write(str(rho))
    f.write(", \"nu\" : ")
    f.write(str(nu))

    f.write(", \"nit\" : ")
    f.write(str(nit))


def results(f,nx,dx,ny,dy,dt,u,v,p,b0):
    print_u (f,",\n \"u_nt\" : [", u, nx, ny, "]")
    print_u (f,",\n \"v_nt\" : [", v, nx, ny, "]")

    print_u (f, ",\n \"b0\" : [", b0, nx, ny, "]")
    b_nt = buildUpB(b0, rho, dt, u, v, dx, dy)

    print_u (f,",\n \"b_nt\" : [", b_nt, nx, ny, "]")
    print_u (f,",\n \"p_nt\" : [", p, nx, ny, "]")
    p_py = presPoisson(p, dx, dy, b_nt)
    print_u (f,",\n \"p_py\" : [", p_py, nx, ny, "]")

    f.write( "}" )



def results1(f,nx,dx,ny,dy,dt,u,v,p,b):
    print_u (f,",\n \"u_nt\" : [", u, nx, ny, "]")
    print_u (f,",\n \"v_nt\" : [", v, nx, ny, "]")
    print_u (f,",\n \"b_nt\" : [", b, nx, ny, "]")
    print_u (f,",\n \"p_nt\" : [", p, nx, ny, "]")
    p_py = presPoisson(p, dx, dy, b)
    print_u (f,",\n \"p_py\" : [", p_py, nx, ny, "]")

    f.write( "}" )



###



def buildUpB(b, rho, dt, u, v, dx, dy):

    b[1:-1,1:-1]=rho*(1/dt*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx)+(v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))-\
                ((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx))**2-\
                2*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dy)*(v[2:,1:-1]-v[0:-2,1:-1])/(2*dx))-\
                ((v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))**2)

    return b


nit=50  ## implicite parameter of presPoisson

def presPoisson(p, dx, dy, b):
    pn = np.empty_like(p)
    pn[:] = p[:]

    for q in range(nit):
                pn[:] = p[:]
                p[1:-1,1:-1] = ((pn[2:,1:-1]+pn[0:-2,1:-1])*dy**2+(pn[1:-1,2:]+pn[1:-1,0:-2])*dx**2)/\
                        (2*(dx**2+dy**2)) -\
                        dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]

                p[-1,:] =p[-2,:]		##dp/dy = 0 at y = 2
                p[0,:] = p[1,:]         ##dp/dy = 0 at y = 0
                p[:,0]=p[:,1]              ##dp/dx = 0 at x = 0
                p[:,-1]=0                      ##p = 0 at x = 2

    return p


def cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))

    for n in range(nt):
        un[:] = u[:]
        vn[:] = v[:]

        b = buildUpB(b, rho, dt, u, v, dx, dy)
        p = presPoisson(p, dx, dy, b)

        u[1:-1,1:-1] = un[1:-1,1:-1]-\
            un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
            vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
            dt/(2*rho*dx)*(p[2:,1:-1]-p[0:-2,1:-1])+\
            nu*(dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+\
            dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2]))

        v[1:-1,1:-1] = vn[1:-1,1:-1]-\
            un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
            vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
            dt/(2*rho*dy)*(p[1:-1,2:]-p[1:-1,0:-2])+\
            nu*(dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])+\
            (dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])))

        u[0,:] = 0
        u[:,0] = 0
        u[:,-1] = 1  ## in last line overwritten below
        v[0,:] = 0
        v[-1,:]=0
        v[:,0] = 0
        v[:,-1] = 0
        u[-1,:] = 0

    return u, v, p, b




nx = 41
ny = 41
nt = 500   ## subsequently overwritten
## nit=50  ## implicite parameter of presPoisson
c = 1      ## unused
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
Y,X = np.meshgrid(y,x)

rho = 1.
nu = .1
dt = .001

# u = np.zeros((ny, nx))
# v = np.zeros((ny, nx))
# p = np.zeros((ny, nx))
# b = np.zeros((ny, nx))

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))
b = np.zeros((ny, nx))
b0 = np.zeros((ny, nx))
nt = 200  ## variable test case parameter

p_py = np.zeros((ny, nx)) ## for testing presPoisson()

### first test case

f = open(''.join(['test-11-', str(nt), '.json']),'w')

print_params(f,nx,dx,ny,dy,nt,dt,rho,nu,nit)

print_u (f, ",\n \"p0\" : [", p, nx, ny, "]")
print_u (f, ",\n \"u0\" : [", u, nx, ny, "]")
print_u (f, ",\n \"v0\" : [", v, nx, ny, "]")

u, v, p, b0 = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)

results(f,nx,dx,ny,dy,dt,u,v,p,b0)

f.close()

fig = plt.figure(figsize=(11,7), dpi=100)
plt.contourf(X,Y,p,alpha=0.5)    ###plotting the pressure field as a contour
plt.colorbar()
plt.contour(X,Y,p)               ###plotting the pressure field outlines
plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
plt.xlabel('X')
plt.ylabel('Y')
fig.show()
#input(" press ENTER to continue ")


### second test case

###  reinitialization

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))
b = np.zeros((ny, nx))
nt = 700  ## variable test case parameter

f = open(''.join(['test-11-', str(nt), '.json']),'w')

print_params(f,nx,dx,ny,dy,nt,dt,rho,nu,nit)

print_u (f, ",\n \"p0\" : [", p, nx, ny, "]")
print_u (f, ",\n \"u0\" : [", u, nx, ny, "]")
print_u (f, ",\n \"v0\" : [", v, nx, ny, "]")


u, v, p, b0 = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)

results(f,nx,dx,ny,dy,dt,u,v,p,b0)

f.close()


fig = plt.figure(figsize=(11,7), dpi=100)
plt.contourf(X,Y,p,alpha=0.5)
plt.colorbar()
plt.contour(X,Y,p)
plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2])
plt.xlabel('X')
plt.ylabel('Y')
fig.show()
#input(" press ENTER to continue ")


### third test case

###  reinitialization

nx = 21
ny = 21

dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
Y,X = np.meshgrid(y,x)

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))
b = np.zeros((ny, nx))
nt = 3  ## variable test case parameter

f = open(''.join(['test-11-', str(nt), '.json']),'w')

print_params(f,nx,dx,ny,dy,nt,dt,rho,nu,nit)

print_u (f, ",\n \"p0\" : [", p, nx, ny, "]")
print_u (f, ",\n \"u0\" : [", u, nx, ny, "]")
print_u (f, ",\n \"v0\" : [", v, nx, ny, "]")


u, v, p, b0 = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)

results(f,nx,dx,ny,dy,dt,u,v,p,b0)

f.close()


fig = plt.figure(figsize=(11,7), dpi=100)
plt.contourf(X,Y,p,alpha=0.5)
plt.colorbar()
plt.contour(X,Y,p)
plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2])
plt.xlabel('X')
plt.ylabel('Y')
fig.show()
input(" press ENTER to continue ")


### fourth test case



def test_cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu, b):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    ## b = np.zeros((ny, nx))

    #for n in range(nt):
    un[:] = u[:]
    vn[:] = v[:]

    ## b = buildUpB(b, rho, dt, u, v, dx, dy)
    ## p = presPoisson(p, dx, dy, b)

    u[1:-1,1:-1] = un[1:-1,1:-1]-\
        un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
        vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
        dt/(2*rho*dx)*(p[2:,1:-1]-p[0:-2,1:-1])+\
        nu*(dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+\
        dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2]))

    v[1:-1,1:-1] = vn[1:-1,1:-1]-\
        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
        dt/(2*rho*dy)*(p[1:-1,2:]-p[1:-1,0:-2])+\
        nu*(dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])+\
        (dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])))

    u[0,:] = 0
    u[:,0] = 0
    u[:,-1] = 1
    v[0,:] = 0
    v[-1,:]=0
    v[:,0] = 0
    v[:,-1] = 0
    u[-1,:] = 0

    return u, v, p, b




###  reinitialization

nx = 21
ny = 21

dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
Y,X = np.meshgrid(y,x)

### take u,v,p,b from third test case

## u = np.zeros((ny, nx))
## v = np.zeros((ny, nx))
## p = np.zeros((ny, nx))
b = np.zeros((ny, nx))  ## and reuse b0 of 3rd test case

nt = 1  ## variable test case parameter / pseudo in this case

f = open(''.join(['test-11-', str(nt), '.json']),'w')

print_params(f,nx,dx,ny,dy,nt,dt,rho,nu,nit)

print_u (f, ",\n \"p0\" : [", p, nx, ny, "]")
print_u (f, ",\n \"u0\" : [", u, nx, ny, "]")
print_u (f, ",\n \"v0\" : [", v, nx, ny, "]")
print_u (f, ",\n \"b0\" : [", b0, nx, ny, "]")

u, v, p, b = test_cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu, b0)
results1(f,nx,dx,ny,dy,dt,u,v,p,b)

f.close()


## fig = plt.figure(figsize=(11,7), dpi=100)
## plt.contourf(X,Y,p,alpha=0.5)
## plt.colorbar()
## plt.contour(X,Y,p)
## plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2])
## plt.xlabel('X')
## plt.ylabel('Y')
## fig.show()
## input(" press ENTER to continue ")
