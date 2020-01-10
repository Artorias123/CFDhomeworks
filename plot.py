#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
filename="out"
rd=open(filename+".dat")
r=rd.readline()
r=[float(a) for a in r.split()]
N=int(r[0])
NS=int(r[1])
rb=r[2]
l=r[3]
#print(N,NS,rb,l)
i=0
val=np.zeros((NS+1,N))
for i in range(0,NS+1):
    r=rd.readline()
    r=[float(a) for a in r.split()]
    val[i]=r
rd.close()
fig, ax = plt.subplots()
xdata, ydata = [], []
pt, =plt.plot([],[],color='black',linestyle='solid',animated=True)
#print(val)
#pt, = plt.plot([], [],'ro-',animated=True)
def init():
    ax.set_xlim(0,N)
    ax.set_ylim(-2,2)
    return pt,
def update(i):
    pt.set_data(range(0,N), val[i])
    return pt,
anim = animation.FuncAnimation(fig, update,range(0,NS+1),interval=100,
init_func=init,blit=True,repeat=True)
anim.save(filename+".gif",writer='imagemagick')
#plt.show()
