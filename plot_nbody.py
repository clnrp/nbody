#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import *

from nbody import *

#------------------------------
def plot_body(x, nb):
    fig = plt.figure()
    ax = Axes3D(fig)
    for n in range(nb):
       xn=[x[i][n] for i in range(int(T/dt))] # posição corpo n 
       xp=[xn[i][0] for i in range(int(T/dt))]
       yp=[xn[i][1] for i in range(int(T/dt))]
       zp=[xn[i][2] for i in range(int(T/dt))]
       ax.plot(xp, yp, zp, zs=0, zdir='z', label='zs=0, zdir=z')
    plt.show()

#------------------------------
def init():
    ax.set_xlim((-D, D))
    ax.set_ylim((-D, D))
    ax.set_zlim((-D, D))
    ax.view_init(30, 0)
    ax.grid(False)
    for linha, ponto in zip(linhas, pontos):
        linha.set_data([], [])
        linha.set_3d_properties([])

        ponto.set_data([], [])
        ponto.set_3d_properties([])
    return linhas+pontos

def animate(i):
    n=0 
    for linha, ponto in zip(linhas, pontos):  
       xn=[xl[j][n] for j in range(i)] # posição corpo n 
       xp=[xn[j][0] for j in range(i)]
       yp=[xn[j][1] for j in range(i)]
       zp=[xn[j][2] for j in range(i)]

       linha.set_data(xp, yp)
       linha.set_3d_properties(zp)

       ponto.set_data(xp[-1:], yp[-1:])
       ponto.set_3d_properties(zp[-1:])
       n+=1

    ax.view_init(30, 0.3 * i)
    fig.canvas.draw()
    return linhas+pontos

#------------------------------
# orbita circular de um satelite
#G=6.674287e-11
#M=5.9722e24
#Tp=10.0
#Dp=((Tp/(2.0*pi))**2*(G*M))**(1.0/3)
#Vp=sqrt(G*M/Dp)
#T=5*Tp
#dt=0.02

# orbita eliptica
t_hora=3600 # 1h = 3600s
G=6.674287e-11 # m^3 * kg^-1 * s^-2
G=G*(t_hora**2) # m^3 * kg^-1 * hora^-2
M=5.9722e24 # kg
e=0.5
P=2 # periodo em horas
K=(4*pi**2)/(G*M)
a=(P**2/K)**(1/3.0)

#eq elipse, 0 no afelio
r = a*(1-e**2)/(1-e*cos(0))

#vis-viva
v = sqrt(G*M*(2/r-1/a))

Dp=r
Vp=v

T=P*3
dt=0.008


D=Dp
print (Dp,Vp)

m=[]
x0=[]
v0=[]

m.append(M) # m1
m.append(5e2) # m2
#m.append(3e2) # m3

x0.append(array([0.0,0.0,0.0])) # x1
x0.append(array([D,0,0])) # x2
#x0.append(npy.array([0.0,D,0])) # x3

v0.append(array([Vp*0.00,Vp*0.00,Vp*0.00])) # v1
v0.append(array([0.0,Vp,0.0])) # v2
#v0.append(npy.array([0.0,0.0,Vp])) # v3

nbody = nbody(G)
[t,x,v,a]=nbody.calculate(m, x0, v0, dt, T, 'rk4')   

dphi=2*pi/250
#r = a*(1-e**2)/(1-e*cos(dphi*125))
#x = r*cos(dphi*125)
#y = r*sin(dphi*125)

sz=len(x)
nd=5
xl=[x[nd*i] for i in range(int(sz/nd))] # reduzir numero de posições
#print x[len(x)-10:],len(x)
print (x[0],x[99],x[100])

fig = plt.figure() #figsize=(10, 10)
ax = fig.add_axes([0.15, 0.15, 0.70, 0.70], projection='3d')

color = plt.cm.jet(linspace(0, 1, len(m)))
linhas = [ax.plot([], [], [], '-', c=c)[0] for c in color]
pontos = [ax.plot([], [], [], 'o', c=c)[0] for c in color]

l=len(xl)
print (l)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=l, interval=50, blit=False, repeat=False)
plt.show()
#anim.save('output.mp4', fps=20)

#vpy=[]
#vpy.append([1,color.red])
#vpy.append([2,color.blue])

#plot_body(x, len(m))
#x.pop(0) remove mais antigo
