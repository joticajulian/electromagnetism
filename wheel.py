from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import json

# Rotation matrices
def R(angle):
    r = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    return r

R60  = R(np.pi/3)
R120 = R(np.pi*2/3)
R180 = R(np.pi)
R240 = R(np.pi*4/3)
R300 = R(np.pi*5/3)

Ic = 1

# rail 1
y1 = np.linspace(-4,0,20)
x1 = 0*y1+3
z1 = 0*y1-1.5
i1 = 0*y1+Ic
c1 = np.vstack((x1,y1,z1)).T
o = np.array([0,0,0])
dl1 = np.vstack((o,np.diff(c1,axis=0)))

# circunference wheel #1
t2a = np.linspace(np.pi*3/2,np.pi*5/2-0.01,20)
t2b = np.linspace(np.pi*3/2-0.01,np.pi/2,20)

y2a = 1.5*np.cos(t2a)
z2a = 1.5*np.sin(t2a)
y2b = 1.5*np.cos(t2b)
z2b = 1.5*np.sin(t2b)

x2 = 0*y2a + 3
i2 = 0*y2a + Ic/2

o = np.array([0,0,0])
c2a = np.vstack((x2,y2a,z2a)).T
c2b = np.vstack((x2,y2b,z2b)).T
dl2a = np.vstack((o,np.diff(c2a,axis=0)))
dl2b = np.vstack((o,np.diff(c2b,axis=0)))

# log curves wheel #1
y3a = np.linspace(1.2593,0,16)
z3a = 0.249616811*np.log(20*y3a+1)
x3 = 0*y3a + 3
i3 = 0*y3a + Ic/6

yz3a = np.vstack((y3a,z3a))
yz3b =  R60.dot(yz3a)
yz3c = R120.dot(yz3a)
yz3d = R180.dot(yz3a)
yz3e = R240.dot(yz3a)
yz3f = R300.dot(yz3a)

c3a = np.vstack((x3,yz3a)).T
c3b = np.vstack((x3,yz3b)).T
c3c = np.vstack((x3,yz3c)).T
c3d = np.vstack((x3,yz3d)).T
c3e = np.vstack((x3,yz3e)).T
c3f = np.vstack((x3,yz3f)).T

o = np.array([0,0,0])
dl3a = np.vstack((o,np.diff(c3a,axis=0)))
dl3b = np.vstack((o,np.diff(c3b,axis=0)))
dl3c = np.vstack((o,np.diff(c3c,axis=0)))
dl3d = np.vstack((o,np.diff(c3d,axis=0)))
dl3e = np.vstack((o,np.diff(c3e,axis=0)))
dl3f = np.vstack((o,np.diff(c3f,axis=0)))

# Axis
x4 = np.linspace(2.9,-2.9,20)
y4 = 0*x4
z4 = 0*x4
i4 = 0*x4+Ic
c4 = np.vstack((x4,y4,z4)).T
o = np.array([0,0,0])+0.1
dl4 = np.vstack((o,np.diff(c1,axis=0)))

# Log curves wheel #2
c5a = np.copy(c3a)[::-1]
c5b = np.copy(c3b)[::-1]
c5c = np.copy(c3c)[::-1]
c5d = np.copy(c3d)[::-1]
c5e = np.copy(c3e)[::-1]
c5f = np.copy(c3f)[::-1]

c5a[:,0] = -3
c5b[:,0] = -3
c5c[:,0] = -3
c5d[:,0] = -3
c5e[:,0] = -3
c5f[:,0] = -3

i5 = i3

o = np.array([0,0,0])
dl5a = np.vstack((o,np.diff(c5a,axis=0)))
dl5b = np.vstack((o,np.diff(c5b,axis=0)))
dl5c = np.vstack((o,np.diff(c5c,axis=0)))
dl5d = np.vstack((o,np.diff(c5d,axis=0)))
dl5e = np.vstack((o,np.diff(c5e,axis=0)))
dl5f = np.vstack((o,np.diff(c5f,axis=0)))

# circunference wheel #2
c6a = np.copy(c2a)[::-1]
c6b = np.copy(c2b)[::-1]

c6a[:,0] = -3
c6b[:,0] = -3

i6 = i2

o = np.array([0,0,0])
dl6a = np.vstack((o,np.diff(c6a,axis=0)))
dl6b = np.vstack((o,np.diff(c6b,axis=0)))

# rail 1
c7 = np.copy(c1)[::-1]
c7[:,0] = -3
i7 = i1
o = np.array([0,0,0])
dl7 = np.vstack((o,np.diff(c7,axis=0)))

c  = np.concatenate((c1,c2a,c2b,c3a,c3b,c3c,c3d,c3e,c3f,c4,c5a,c5b,c5c,c5d,c5e,c5f,c6a,c6b,c7))
dl  = np.concatenate((dl1,dl2a,dl2b,dl3a,dl3b,dl3c,dl3d,dl3e,dl3f,dl4,dl5a,dl5b,dl5c,dl5d,dl5e,dl5f,dl6a,dl6b,dl7))
I  = np.concatenate((i1,i2,i2,i3,i3,i3,i3,i3,i3,i4,i5,i5,i5,i5,i5,i5,i6,i6,i7))
I0 = np.concatenate((0*i1,0*i2,0*i2,i3,i3,i3,i3,i3,i3,0*i4,i5,i5,i5,i5,i5,i5,0*i6,0*i6,0*i7))

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
    
f = F(c,dl,I0)
#jsondata = {'circuit':c.tolist(),'force':f.tolist()}
#with open("datawheel.json", 'w') as outfile:
#   json.dump(jsondata, outfile)

#with open("datawheel.json", 'r') as file:
#    jsondata = json.load(file) 

f[I0==0] = 0

# 3d figure
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(c1[:,0],c1[:,1],c1[:,2],color='r')
ax.plot(c2a[:,0],c2a[:,1],c2a[:,2],color='r')
ax.plot(c2b[:,0],c2b[:,1],c2b[:,2],color='r')
ax.plot(c3a[:,0],c3a[:,1],c3a[:,2],color='r')
ax.plot(c3b[:,0],c3b[:,1],c3b[:,2],color='r')
ax.plot(c3c[:,0],c3c[:,1],c3c[:,2],color='r')
ax.plot(c3d[:,0],c3d[:,1],c3d[:,2],color='r')
ax.plot(c3e[:,0],c3e[:,1],c3e[:,2],color='r')
ax.plot(c3f[:,0],c3f[:,1],c3f[:,2],color='r')
ax.plot(c4[:,0],c4[:,1],c4[:,2],color='r')
ax.plot(c5a[:,0],c5a[:,1],c5a[:,2],color='r')
ax.plot(c5b[:,0],c5b[:,1],c5b[:,2],color='r')
ax.plot(c5c[:,0],c5c[:,1],c5c[:,2],color='r')
ax.plot(c5d[:,0],c5d[:,1],c5d[:,2],color='r')
ax.plot(c5e[:,0],c5e[:,1],c5e[:,2],color='r')
ax.plot(c5f[:,0],c5f[:,1],c5f[:,2],color='r')
ax.plot(c6a[:,0],c6a[:,1],c6a[:,2],color='r')
ax.plot(c6b[:,0],c6b[:,1],c6b[:,2],color='r')
ax.plot(c7[:,0],c7[:,1],c7[:,2],color='r')

subgroup = range(0,c.shape[0],3)
ch = c[subgroup]
fh = f[subgroup]
ax.quiver(ch[:,0],ch[:,1],ch[:,2],fh[:,0],fh[:,1],fh[:,2],length=0.4,normalize=True)
plt.title('Forces')

plt.show()
