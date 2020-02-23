#!/usr/bin/python
# -*- coding: UTF-8 -*-

#from pylab import *
#from scitools.easyviz.gnuplot_ import *
import pickle
import sys
import time
from random import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import *
import matplotlib.animation as animation

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
    ax.set_zlim((-D+8, D-8))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(30, 120)
    #ax.grid(False)
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

    #ax.view_init(30, 0.3 * i)
    fig.canvas.draw()
    return linhas+pontos

def val(ref): # gera valores aleatórios
    return randrange(-100,100)/100.*ref

def mult_stars(Pp,Vp,Mp,Dp,N,sn):
    [m,x0,v0]=sn
    # corpo M0
    m.append(Mp*.8) 
    x0.append(Pp) 
    v0.append(Vp) 

    # outros corpos
    for i in range(1,N): 
        m.append(0.15*Mp*abs(val(1))) # m
        x0_a=array([val(Dp),val(Dp),val(Dp)])+Pp
        x0.append(x0_a) # posição

        rp=x0_a-Pp # posição em relação a M0
        R=linalg.norm(rp)
        Vn=sqrt(G*Mp/R)
        
        # o vetor velocidade é ortogonal a rp
        vt=array([rp[1], -rp[0], 0])
        vtn=linalg.norm(vt)
        vel=Vn*vt/vtn+0.0*Vp
        v0.append(vel) # velocidade
    return [m,x0,v0]

def mape(x1, x2, N, ini):
    err=0.0
    errx=erry=errz=0
    for i in range(ini,N): 
        x1p=[x1[j][i][0] for j in range(len(x1))]
        y1p=[x1[j][i][1] for j in range(len(x1))]
        z1p=[x1[j][i][2] for j in range(len(x1))]
    
        x2p=[x2[j][i][0] for j in range(len(x2))]
        y2p=[x2[j][i][1] for j in range(len(x2))]
        z2p=[x2[j][i][2] for j in range(len(x2))]

        for j in range(len(x2p)):
            if(abs(x1p[j])<1):x1p[j]=float(1.0)
            if(abs(y1p[j])<1):y1p[j]=float(1.0)
            if(abs(z1p[j])<1):z1p[j]=float(1.0)
            if(abs(x2p[j])<1):x2p[j]=float(1.0)
            if(abs(y2p[j])<1):y2p[j]=float(1.0)
            if(abs(z2p[j])<1):z2p[j]=float(1.0)

        errx=[abs((x2p[j]-x1p[j])/x1p[j]*100) for j in range(1,len(x1p))]
        erry=[abs((y2p[j]-y1p[j])/y1p[j]*100) for j in range(1,len(y1p))]
        errz=[abs((z2p[j]-z1p[j])/z1p[j]*100) for j in range(1,len(z1p))]
        err+=(sum(errx)+sum(erry)+sum(errz))/(3*len(errx))
    err=err/(N-ini)
    return err

#------------------------------
#if __name__ == "__main__":

larg=len(sys.argv);
cmd=''
nfile=''
if(larg>1):
   cmd=sys.argv[1]
   nfile=sys.argv[2]
   print(larg,cmd,nfile)

#with open('', 'wb') as f:
#    pickle.dump(mlist, f)

d_ua=149597870700
t_ano=31536000
m_sol=1.989e30
G=6.674287e-11  # m^3 * kg^-1 * s^-2
G=G*(t_ano**2)*(m_sol)/(d_ua**3) # ua^3 * msol^-1 * ano^-2
Mp=3 # 3 massa do sol
N=5
T=60 # 30 anos
dt=0.005
Dp=12 # 4 ua 
D=Dp*3

m=[]
x0=[]
v0=[]
sn=[m,x0,v0]
mlist=[]

if(larg==3 and cmd=='r'):
   with open(nfile, 'rb') as f:
      mlist = pickle.load(f)
      m=mlist[0]
      x0=mlist[1]
      v0=mlist[2] 
      print(mlist)
else:
   Pp=array([0,0,0])
   Vp=array([0.8,0,0])
   [m,x0,v0]=mult_stars(Pp,Vp,Mp,Dp,N,sn)
   mlist=[m,x0,v0]
   if(larg==3 and cmd=='w'):
      print('gravar!')
      with open(nfile, 'wb') as f:
         pickle.dump(mlist, f)
   else:
      with open('nb_m.l', 'wb') as f:
         pickle.dump(mlist, f)

print(m,x0,v0)

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

err1=mape(x2, x1, len(m), 0)
err2=mape(x2, x1, len(m), 0)
print('err1:'+str(err1))
print('err2:'+str(err2))

sz=len(x1)
nd=25
xl=[x2[nd*i] for i in range(int(sz/nd))] # reduzir numero de posições
#print x[len(x)-10:],len(x)
print(x1[0],x1[99],x1[100])

fig = plt.figure() #figsize=(10, 10)
ax = fig.add_axes([0.15, 0.15, 0.70, 0.70], projection='3d')

color = plt.cm.jet(linspace(0, 1, len(m)))
linhas = [ax.plot([], [], [], '-', c=c)[0] for c in color]
pontos = [ax.plot([], [], [], 'o', c=c)[0] for c in color]

l=len(xl)
print(t1x,t2x)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=l, interval=50, blit=False, repeat=False)
#anim.save('ani_nb_mult.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()

#vpy=[]
#vpy.append([1,color.red])
#vpy.append([2,color.blue])

#plot_body(x, len(m))
#x.pop(0) remove mais antigo
