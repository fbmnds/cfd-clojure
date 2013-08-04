import numpy as np                 #loading our favorite library
import matplotlib.pyplot as plt    #and the useful plotting library
from IPython.core.display import clear_output #used for inline animation

nx = 41
dx = 2./(nx-1)
nt = 20    #the number of timesteps we want to calculate
nu = 0.3   #the value of viscosity
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = sigma*dx**2/nu #dt is defined using sigma ... more later!


u = np.ones(nx)      #a numpy array with nx elements all equal to 1.
u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s

un = np.ones(nx) #our placeholder array, un, to advance the solution in time

for n in range(nt):  #iterate through time
    un[:] = u[:] ##copy the existing values of u into un
    for i in range(1,nx-1):
        u[i] = un[i] + nu*dt/dx**2*(un[i+1]-2*un[i]+un[i-1])

print (u)
plt.plot(np.linspace(0,2,nx), u)
plt.show()


# # Remember: comments in python are denoted by the pound sign
# import numpy as np                 #here we load numpy, calling it 'np' from now on
# import matplotlib.pyplot as plt    #here we load matplotlib, calling it 'plt'
# import time, sys                   #and load some utilities
# from IPython.core.display import clear_output #used for inline animation
#
# nx = 41  # try changing this number from 41 to 81 and Run All ... what happens?
# dx = 2./(nx-1)
# nt = 20    #nt is the number of timesteps we want to calculate
# dt = .025  #dt is the amount of time each timestep covers (delta t)
# c = 1.      #assume wavespeed of c = 1
#
# u = np.ones(nx)      #numpy function ones()
# u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s
# print (u)
#
# plt.plot(np.linspace(0,2,nx), u)
# plt.show()
#
# un = np.ones(nx) #initialize a temporary array
#
# for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
#     un[:] = u[:] ##copy the existing values of u into un
#     #for i in range(1,nx): ## you can try commenting this line and...
#     for i in range(nx): ## ... uncommenting this line and see what happens!
#         u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1])
#
# print (u)
#
# plt.plot(np.linspace(0,2,nx),u)
# plt.show()
#
#
# nx = 81  # try changing this number from 41 to 81 and Run All ... what happens?
# dx = 2./(nx-1)
# nt = 20    #nt is the number of timesteps we want to calculate
# dt = .001  #dt is the amount of time each timestep covers (delta t)
# c = 1.      #assume wavespeed of c = 1
#
# u = np.ones(nx)      #numpy function ones()
# u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s
# print (u)
#
# plt.plot(np.linspace(0,2,nx), u)
# plt.show()
#
# un = np.ones(nx) #initialize a temporary array
#
# for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
#     un[:] = u[:] ##copy the existing values of u into un
#     #for i in range(1,nx): ## you can try commenting this line and...
#     for i in range(nx): ## ... uncommenting this line and see what happens!
#         u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1])
#
# print (u)
#
#
# plt.plot(np.linspace(0,2,nx),u)
# plt.show()
