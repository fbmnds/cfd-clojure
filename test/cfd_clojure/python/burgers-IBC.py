import numpy as np
import matplotlib.pyplot as plt

###variable declarations
nx = 101
nt = 100
dx = 2*np.pi/(nx-1)
nu = .07
dt = dx*nu
un = np.empty(nx)

x = np.linspace(0, 2*np.pi, nx)
print ('nx, nt, dx, nu, dt: ', nx, nt, dx, nu, dt)
print ('\nx : ')
print (x)

def exp(x):
    return np.exp(x)

def burgers_u0 (t,x,nu):
    return -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 12.5663706143592)*exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1)))) + 4

u = burgers_u0 (0,x,nu)
print ('\nBurgers u0 : ')
print (u)


for n in range(nt):
    un[:]=u[:]
    for i in range(nx-1):
        u[i] = un[i] - un[i] * dt/dx *(un[i] - un[i-1]) + nu*dt/dx**2*\
                (un[i+1]-2*un[i]+un[i-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
                (un[0]-2*un[-1]+un[-2])

u_analytical = np.asarray([burgers_u0(nt*dt, xi, nu) for xi in x])

print ('\nu : ')
print (u)
print ('\nu(analytical) : ')
print (u_analytical)



nt = 1
u = burgers_u0 (0,x,nu)
for n in range(nt):
    un[:]=u[:]
    for i in range(nx-1):
        u[i] = un[i] - un[i] * dt/dx *(un[i] - un[i-1]) + nu*dt/dx**2*\
                (un[i+1]-2*un[i]+un[i-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
                (un[0]-2*un[-1]+un[-2])
    print('\ un[0] un[-1]', un[0], un[-1])

print ('\nnt : ', nt)
print ('\nu : ')
print (u)


nt = 2
u = burgers_u0 (0,x,nu)
for n in range(nt):
    un[:]=u[:]
    for i in range(nx-1):
        u[i] = un[i] - un[i] * dt/dx *(un[i] - un[i-1]) + nu*dt/dx**2*\
                (un[i+1]-2*un[i]+un[i-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
                (un[0]-2*un[-1]+un[-2])
    print('\ un[0] un[-1]', un[0], un[-1])

print ('\nnt : ', nt)
print ('\nu : ')
print (u)
