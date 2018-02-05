from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import json

# spinning coil - geometry

# vertical line
z1 = np.linspace(5,2.1,29)
x1 = 0*z1
y1 = 0*z1

# axial line
x2 = np.linspace(0.1,1.4,15)
y2 = 0*x2
z2 = 0*x2 + 2

# coil
a = np.linspace(0,2*np.pi*5,471)
x3 = 1.5*np.cos(a)
y3 = 1.5*np.sin(a)
z3 = np.linspace(2,0,471)

cx = np.concatenate((x1,x2,x3))
cy = np.concatenate((y1,y2,y3))
cz = np.concatenate((z1,z2,z3))

c = np.vstack((cx,cy,cz)).T

def B(circuit,I,x,y,z):
    bx = 0*x
    by = 0*y
    bz = 0*z
    
    mu = 1.26 * 10**(-6)
    last_p = circuit[0]    
    for p in circuit:
        dl = p - last_p
        last_p = p
        for i in range(0,x.shape[0]):
            for j in range(0,x.shape[1]):
                for k in range(0,x.shape[2]):
                    r = np.array([x[j,i,k]-p[0] , y[j,i,k]-p[1] , z[j,i,k]-p[2]])                    
                    norm_r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2);
                    if norm_r<=0.01:
                        continue                    
                    dB = np.cross(dl,r)/(norm_r**3)
                    bx[j,i,k] += dB[0]
                    by[j,i,k] += dB[1]
                    bz[j,i,k] += dB[2]        
    bx *= mu*I/(4*np.pi)
    by *= mu*I/(4*np.pi)
    bz *= mu*I/(4*np.pi)
    return bx,by,bz
    
def F(circuit,I):
    force = 0*circuit    
    mu = 1.26 * 10**(-6)
    last_p = circuit[0]
    for p in circuit:
        dl = p - last_p
        last_p = p
        for i in range(0,circuit.shape[0]):
            r = circuit[i]-p
            norm_r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2);
            if norm_r<=0.01:
                continue                    
            dB = np.cross(dl,r)/(norm_r**3)
            dl2 = np.array([0,0,0])
            if i>0:
                dl2 = circuit[i]-circuit[i-1]
            dF = np.cross(dl2,dB)
            force[i] += dF
    force *= mu*I*I/(4*np.pi)
    return force

# x
x = np.linspace(-1.8,1.8,8)
y = np.linspace(-1.8,1.8,8)
z = np.linspace(0,5,6)
x,y,z = np.meshgrid(x,y,z)
bx,by,bz = B(c,1,x,y,z)
f = F(c,1)

#jsondata = {'circuit':c.tolist(),'force':f.tolist(),'magfield':{'x':bx.tolist(),'y':by.tolist(),'z':bz.tolist()},'space':{'x':x.tolist(),'y':y.tolist(),'z':z.tolist()}}
# with open("dataCoil.json", 'w') as outfile:
#    json.dump(jsondata, outfile)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(x,y,z,bx,by,bz,color='b',length=0.5,normalize=True)
ax.plot(c[:,0],c[:,1],c[:,2],color='r')
plt.title('Magnetic Field')

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(c[:,0],c[:,1],c[:,2],color='r')

subgroup = range(0,c.shape[0],5)
c = c[subgroup]
f = f[subgroup]
ax.quiver(c[:,0],c[:,1],c[:,2],f[:,0],f[:,1],f[:,2],length=0.5,normalize=True)
plt.title('Forces')

plt.show()
        


