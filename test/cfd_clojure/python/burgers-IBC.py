import numpy as np

###variable declarations
nx = 101
nt = 100
dx = 2*np.pi/(nx-1)
nu = .07
dt = dx*nu

x = np.linspace(0, 2*np.pi, nx)
print (x)

def exp(x):
    return np.exp(x)

def burgers_u0 (t,x,nu):
    return -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 12.5663706143592)*exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1)))) + 4

print (burgers_u0 (0,x,nu))
