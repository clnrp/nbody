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
    global R
    ax.set_xlim((-D*.2, D))
    ax.set_ylim((-D*.2, D*1.3))
    ax.set_zlim((R, R+h*1.3))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(30, 0)
    #ax.grid(False)
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

    #ax.view_init(30, 0.3*i)
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
#t_hora=3600 # 1h = 3600s
#G=6.674287e-11 # m^3 * kg^-1 * s^-2
#G=G*(t_hora**2) # m^3 * kg^-1 * hora^-2
#M=5.9722e24 # kg
#e=0.5
#P=10 # periodo em horas
#K=(4*pi**2)/(G*M)
#a=(P**2/K)**(1/3.0)

#eq elipse, 0 no afelio
#r = a*(1-e**2)/(1-e*cos(0))

#vis-viva
#v = sqrt(G*M*(2/r-1/a))

G=6.674287e-11 # m^3 * kg^-1 * s^-2
M=5.9722e24 # kg
R=6400e3
V=50 # m/s
phi=pi/3 # ang rad
g=G*M/R**2
T=2*V*sin(phi)/g # tempo (s)
dt=0.01
D=V*cos(phi)*T
h=(V*sin(phi))**2/(2*g) # altura em metros
Vy=V*cos(phi)
Vz=V*sin(phi)

m=[]
x0=[]
v0=[]

m.append(M) # m1
m.append(5e2) # m2

x0.append(array([0.0,0.0,0.0])) # x1
x0.append(array([0,0,R])) # x2

v0.append(array([0.0,0.0,0.0])) # v1
v0.append(array([0.0,Vy,Vz])) # v2

nbody = nbody(G)
tm1 = time.time()
[t1,x1,v1,a1]=nbody.calculate(m, x0, v0, dt, T, 'euler')  
tm2 = time.time()
[t2,x2,v2,a2]=nbody.calculate(m, x0, v0, dt, T, 'rk4')   
tm3 = time.time()
t1x=tm2-tm1
t2x=tm3-tm2

#[t1,x1,v1,a1]=nbody.calculate(m, x0, v0, dt, T, 'euler')  
#[t2,x2,v2,a2]=nbody.calculate(m, x0, v0, dt, T, 'rk4')   


y1p=[x1[i][1][1] for i in range(len(x1))]
z1p=[x1[i][1][2] for i in range(len(x1))]
h=max(z1p)
ip=z1p.index(h)
tp=t1[ip]
h=h-R
d=2*tp*V*cos(phi)
print(ip,h,tp,d,g)
#t=v0z/g=v*sin(phi)/g

x3=map(list.__add__, x1, x2)
sz=len(x3)
nd=1
xl=[x3[nd*i] for i in range(int(sz/nd))] # reduzir numero de posições

y2p=[x2[i][1][1] for i in range(len(x2))]
z2p=[x2[i][1][2] for i in range(len(x2))]

ye=[Vy*t1[i] for i in range(len(t1))]
ze=[6400000+Vz*t1[i]-.5*g*t1[i]**2 for i in range(len(t1))]

erry1=[(y1p[i]-ye[i])/ye[i]*100 for i in range(1,len(ye))]
errz1=[(z1p[i]-ze[i])/ze[i]*100 for i in range(1,len(ze))]
err1=(sum(erry1)/len(erry1)+sum(errz1)/len(errz1))/2

erry2=[(y2p[i]-ye[i])/ye[i]*100 for i in range(1,len(ye))]
errz2=[(z2p[i]-ze[i])/ze[i]*100 for i in range(1,len(ze))]
err2=(sum(erry2)/len(erry2)+sum(errz2)/len(errz2))/2

print(z1p[len(y1p)-5:],z2p[len(y2p)-5:],ze[len(ye)-5:])
print('err1:'+str(err1))
print('err2:'+str(err2))

fig = plt.figure() #figsize=(10, 10)
ax = fig.add_axes([0.15, 0.15, 0.70, 0.70], projection='3d')

color = plt.cm.jet(linspace(0, 1, 2*len(m)))
linhas = [ax.plot([], [], [], '-', c=c)[0] for c in color]
pontos = [ax.plot([], [], [], 'o', c=c)[0] for c in color]

l=len(xl)
print(t1x,t2x)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=l, interval=50, blit=False, repeat=False)
plt.show()

#anim.save('output.mp4', fps=20)

#vpy=[]
#vpy.append([1,color.red])
#vpy.append([2,color.blue])

#plot_body(x, len(m))
#x.pop(0) remove mais antigo
