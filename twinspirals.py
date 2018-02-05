from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import json

# Rotation matrices
R90  = np.array([[ 0,-1],[ 1, 0]])
R180 = np.array([[-1, 0],[ 0,-1]])
R270 = np.array([[ 0, 1],[-1, 0]])

Ic = 1

# vertical line
z1 = np.linspace(5,2.8,20)
x1 = 0*z1
y1 = 0*z1

c1 = np.vstack((x1,y1,z1)).T
o = np.array([0,0,0])
dl1 = np.vstack((o,np.diff(c1,axis=0)))
i1 = 0*z1+Ic

# log curves top
x2a = np.linspace(0,1.2593,16)
y2a = 0.249616811*np.log(20*x2a+1)
z2 = 0*x2a + 2.7

i2 = 0*x2a + Ic/4

xy2a = np.vstack((x2a,y2a))
xy2b =  R90.dot(xy2a)
xy2c = R180.dot(xy2a)
xy2d = R270.dot(xy2a)

c2a = np.vstack((xy2a,z2)).T
c2b = np.vstack((xy2b,z2)).T
c2c = np.vstack((xy2c,z2)).T
c2d = np.vstack((xy2d,z2)).T

o = c2a[0]-c1[-1]
dl2a = np.vstack((o,np.diff(c2a,axis=0)))
dl2b = np.vstack((o,np.diff(c2b,axis=0)))
dl2c = np.vstack((o,np.diff(c2c,axis=0)))
dl2d = np.vstack((o,np.diff(c2d,axis=0)))


# log curves bottom
RA = np.array([[0.9123,0.4095],[-0.4095,0.9123]])
yx2a = np.vstack((y2a[::-1],x2a[::-1]))

xy3a = RA.dot(yx2a)
xy3b = R90.dot(xy3a)
xy3c = R180.dot(xy3a)
xy3d = R270.dot(xy3a)

z3 = 0*x2a + 2.3
i3 = i2

c3a = np.vstack((xy3a,z3)).T
c3b = np.vstack((xy3b,z3)).T
c3c = np.vstack((xy3c,z3)).T
c3d = np.vstack((xy3d,z3)).T

oa = c3a[0]-c2a[-1]
ob = c3b[0]-c2b[-1]
oc = c3c[0]-c2c[-1]
od = c3d[0]-c2d[-1]
dl3a = np.vstack((oa,np.diff(c3a,axis=0)))
dl3b = np.vstack((ob,np.diff(c3b,axis=0)))
dl3c = np.vstack((oc,np.diff(c3c,axis=0)))
dl3d = np.vstack((od,np.diff(c3d,axis=0)))

# vertical line
z4 = np.linspace(2.2,0,20)
x4 = 0*z1
y4 = 0*z1

c4 = np.vstack((x4,y4,z4)).T
o = c4[0]-c3a[-1]
dl4 = np.vstack((o,np.diff(c4,axis=0)))
i4 = 0*z4 + Ic

c  = np.concatenate((c1,c2a,c2b,c2c,c2d,c3a,c3b,c3c,c3d,c4))
dl = np.concatenate((dl1,dl2a,dl2b,dl2c,dl2d,dl3a,dl3b,dl3c,dl3d,dl4))
I  = np.concatenate((i1,i2,i2,i2,i2,i3,i3,i3,i3,i4))

def B(circuit,dl,I,x,y,z):
    bx = 0*x
    by = 0*y
    bz = 0*z
    
    mu = 1.26 * 10**(-6)    
    for m in range(0,circuit.shape[0]):
        for i in range(0,x.shape[0]):
            for j in range(0,x.shape[1]):
                for k in range(0,x.shape[2]):
                    P = np.array([x[j,i,k] , y[j,i,k] , z[j,i,k]])
                    r = P - circuit[m]
                    norm_r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2);
                    if norm_r<=0.01:
                        continue                    
                    dB = np.cross(dl[m],r)/(norm_r**3)                    
                    bx[j,i,k] += dB[0]*I[m]/(4*np.pi)
                    by[j,i,k] += dB[1]*I[m]/(4*np.pi)
                    bz[j,i,k] += dB[2]*I[m]/(4*np.pi)
    return bx,by,bz
 
 
def F(circuit,dl,I):
    force = 0*circuit    
    mu = 1.26 * 10**(-6)
    for m in range(0,circuit.shape[0]):
        for i in range(0,circuit.shape[0]):
            r = circuit[i]-circuit[m]
            norm_r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2);
            if norm_r<=0.01:
                continue                    
            dB = np.cross(dl[m],r)/(norm_r**3)
            dF = np.cross(dl[i],dB)
            force[i] += dF*mu*I[m]*I[i]/(4*np.pi)    
    return force
    
x = np.linspace(-1.8,1.8,8)
y = np.linspace(-1.8,1.8,8)
z = np.linspace(0,5,6)
x,y,z = np.meshgrid(x,y,z)
bx,by,bz = B(c,dl,I,x,y,z)
f = F(c,dl,I)

#jsondata = {'circuit':c.tolist(),'force':f.tolist(),'magfield':{'x':bx.tolist(),'y':by.tolist(),'z':bz.tolist()},'space':{'x':x.tolist(),'y':y.tolist(),'z':z.tolist()}}

#jsondata['circuit-details'] = {'c1':c1.tolist(),'c2':c2.tolist(),'cb1':cb1.tolist(),'cb2':cb2.tolist(),'cc1':cc1.tolist(),'cc2':cc2.tolist(),'cd1':cd1.tolist(),'cd2':cd2.tolist()}

#with open("datatwinspirals.json", 'w') as outfile:
#   json.dump(jsondata, outfile) 
   
   
#with open("datatwinspirals.json", 'r') as file:
#    jsondata = json.load(file)
    
#c = np.array(jsondata['circuit'])
#f = 6e7*np.array(jsondata['force'])
#x = np.array(jsondata['space']['x'])
#y = np.array(jsondata['space']['y'])
#z = np.array(jsondata['space']['z'])
#bx = np.array(jsondata['magfield']['x'])
#by = np.array(jsondata['magfield']['y'])
#bz = np.array(jsondata['magfield']['z'])
f = 6e7*f

# 3d figure
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(x,y,z,bx,by,bz,color='b',length=1,normalize=True)
ax.plot(c1[:,0],c1[:,1],c1[:,2],color='r')
ax.plot(c2a[:,0],c2a[:,1],c2a[:,2],color='r')
ax.plot(c2b[:,0],c2b[:,1],c2b[:,2],color='r')
ax.plot(c2c[:,0],c2c[:,1],c2c[:,2],color='r')
ax.plot(c2d[:,0],c2d[:,1],c2d[:,2],color='r')
ax.plot(c3a[:,0],c3a[:,1],c3a[:,2],color='orange')
ax.plot(c3b[:,0],c3b[:,1],c3b[:,2],color='orange')
ax.plot(c3c[:,0],c3c[:,1],c3c[:,2],color='orange')
ax.plot(c3d[:,0],c3d[:,1],c3d[:,2],color='orange')
ax.plot(c4[:,0],c4[:,1],c4[:,2],color='orange')
plt.title('Magnetic Field')

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(c1[:,0],c1[:,1],c1[:,2],color='r')
ax.plot(c2a[:,0],c2a[:,1],c2a[:,2],color='r')
ax.plot(c2b[:,0],c2b[:,1],c2b[:,2],color='r')
ax.plot(c2c[:,0],c2c[:,1],c2c[:,2],color='r')
ax.plot(c2d[:,0],c2d[:,1],c2d[:,2],color='r')
ax.plot(c3a[:,0],c3a[:,1],c3a[:,2],color='orange')
ax.plot(c3b[:,0],c3b[:,1],c3b[:,2],color='orange')
ax.plot(c3c[:,0],c3c[:,1],c3c[:,2],color='orange')
ax.plot(c3d[:,0],c3d[:,1],c3d[:,2],color='orange')
ax.plot(c4[:,0],c4[:,1],c4[:,2],color='orange')
ax.quiver(c[:,0],c[:,1],c[:,2],f[:,0],f[:,1],f[:,2])
plt.title('Forces')

plt.show()
