import numpy as np

nx = 4
ny = 5

print ("nx, ny : ", nx, ny)


u = np.ones((ny,nx))
for i in range(ny):
    for j in range(nx):
        u[i,j] = i + j/10.

print ("u ")
print (u)

print ("u[1:,1:] ; A :except-rows first :except-cols first ")
print (u[1:,1:])

print ("u[0:-1,1:] ; B :except-rows last :except-cols first ")
print (u[0:-1,1:])

print ("u[1:,0:-1] ; C :except-rows first :except-cols last ")
print (u[1:,0:-1])

print ("u[1:-1,1:-1] ; E :except-rows first  :except-rows last :except-cols first :except-cols last ")
print (u[1:-1,1:-1])

print ("u[2:,1:-1] ; F :except-rows first  :except-rows second :except-cols first :except-cols last ")
print (u[2:,1:-1])

print ("u[0:-2,1:-1] ; G :except-rows prev-last  :except-rows last :except-cols first :except-cols last  ")
print (u[0:-2,1:-1])

print ("u[1:-1,2:] ; H :except-rows first  :except-rows last :except-cols first :except-cols second   ")
print (u[1:-1,2:])

print ("u[1:-1,0:-2] ; J :except-rows first  :except-rows last :except-cols prev-last :except-cols last   ")
print (u[1:-1,0:-2])

print("u[0,:] ; first row ")
print(u[0,:])

print("u[-1,:] ; last row ")
print(u[-1,:])

print("u[:,0] ; first col ")
print(u[:,0])

print("u[:,-1] ; last col ")
print(u[:,-1])


print ("u**2 ")
print (u**2)
