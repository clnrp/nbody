#!/usr/bin/python
# -*- coding: UTF-8 -*-

from random import *
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
    pmax=100
    if(i>pmax): 
       p=pmax 
    else: 
       p=i
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

def val(ref): # gera valores aleatórios
    return randrange(-100,100)/100.*ref

def mult_stars(Pp,Vp,Mp,Dp,N):
    global m,x0,v0
    # corpo M0
    m=[]
    m.append(Mp) 
    x0.append(Pp) 
    v0.append(Vp) 

    # outros corpos
    for i in range(1,N): 
        m.append(0.1*Mp*abs(val(1))) # m
        x0_a=array([val(Dp),val(Dp),val(Dp)])+Pp
        x0.append(x0_a) # posição

        rp=x0_a-Pp # posição em relação a M0
        R=linalg.norm(rp)
        Vn=sqrt(G*Mp/R)
        
        # o vetor velocidade é ortogonal a rp
        vt=array([rp[1], -rp[0], 0])
        vtn=linalg.norm(vt)
        vel=Vn*vt/vtn+0.7*Vp
        v0.append(vel) # velocidade

#------------------------------
d_ua=149597870700
t_ano=31536000
m_sol=1.989e30
G=6.674287e-11  # m^3 * kg^-1 * s^-2
G=G*(t_ano**2)*(m_sol)/(d_ua**3) # ua^3 * msol^-1 * ano^-2
Mp=1 # 1 massa do sol
N1=5
N2=6
T=30 # 30 anos
dt=0.02
Dp=4 # 4 ua 
D=Dp*3

m=[]
x0=[]
v0=[]

# sistema 1
Pp=array([-20,0,0])
Vp=array([0.8,0,0])
mult_stars(Pp,Vp,Mp,Dp,N1)

# sistema 2
Pp=array([20,0,0])
Vp=array([-0.8,0,0])
mult_stars(Pp,Vp,Mp,Dp,N2)

for i in range(len(m)):
    print('m:'+str(m[i]),'x0:'+str(x0[i]),'v0:'+str(v0[i]))

nbody = nbody(G)
[t,x,v,a]=nbody.calculate(m, x0, v0, dt, T, 'euler')   

sz=len(x)
nd=5
xl=[x[nd*i] for i in range(sz/nd)] # reduzir numero de posições
#print x[len(x)-10:],len(x)
print(x[0],x[99],x[100])

fig = plt.figure() #figsize=(10, 10)
ax = fig.add_axes([0.15, 0.15, 0.70, 0.70], projection='3d')

color = plt.cm.jet(linspace(0, 1, len(m)))
linhas = [ax.plot([], [], [], '-', c=c)[0] for c in color]
pontos = [ax.plot([], [], [], 'o', c=c)[0] for c in color]

l=len(xl)
print(l)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=l, interval=50, blit=False, repeat=False)
plt.show()
#anim.save('output.mp4', fps=20)

#vpy=[]
#vpy.append([1,color.red])
#vpy.append([2,color.blue])

#plot_body(x, len(m))
#x.pop(0) remove mais antigo
