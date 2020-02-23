#!/usr/bin/python
# -*- coding: UTF-8 -*-
from numpy import *

class nbody(object):  

    # class constructor
    def __init__(self, G=6.674287e-11): 
        self.G = G 
        print ('G=' + str(G))
    
    # calculate acceleration
    def f(self, i, m, x): 
        G = self.G
        ax = ay = az = 0
        for j in range(self.nb):
           if(i != j):  # or multiplying by (1-Delta[i][j]) where Delta is the Kronecker delta function
              ax += -G * m[j] * (x[i][0] - x[j][0]) / (((x[i][0] - x[j][0]) ** 2 + (x[i][1] - x[j][1]) ** 2 + (x[i][2] - x[j][2]) ** 2) ** (3.0 / 2))
              ay += -G * m[j] * (x[i][1] - x[j][1]) / (((x[i][0] - x[j][0]) ** 2 + (x[i][1] - x[j][1]) ** 2 + (x[i][2] - x[j][2]) ** 2) ** (3.0 / 2))
              az += -G * m[j] * (x[i][2] - x[j][2]) / (((x[i][0] - x[j][0]) ** 2 + (x[i][1] - x[j][1]) ** 2 + (x[i][2] - x[j][2]) ** 2) ** (3.0 / 2))
        a = array([ax, ay, az])
        return a

    # calculates speed and position using Euler's method
    def euler(self, f, m, x0, v0, dt): 
        v = list(self.vet0_l)
        x = list(self.vet0_l)
        a = list(self.vet0_l)
        for i in range(self.nb):
           a[i] = f(i, m, x0)
           x[i] = x0[i] + dt * v0[i]  # x=x0+v*dt
           v[i] = v0[i] + dt * a[i]  # v=v0+a*dt
        return [x, v, a]
    
    # calculates speed and position using improved Euler method
    def eulerm(self, f, m, x0, v0, dt): 
        v = list(self.vet0_l)
        x = list(self.vet0_l)
        a = list(self.vet0_l)
        for i in range(self.nb):
           x1 = list(x0)
           a1 = f(i, m, x0)  # a1=f(m,x0)
           x1[i] = x0[i] + dt * v0[i]  # x1=x0+v0*dt
           v1 = v0[i] + dt * a1  # v1=v0+a1*dt
 
           a[i] = 0.5 * (a1 + f(i, m, x1))  # a=1/2*(f(m,x0)+f(m,x1))
           x[i] = x0[i] + dt * (0.5 * (v0[i] + v1))  # x=x0+1/2*(v0+v1)*dt
           v[i] = v0[i] + dt * a[i]  # v=v0+a*dt
        return [x, v, a]
    
    # calculates speed and position using Runge Kutta method of 4th order
    def rk4(self, f, m, x0, v0, dt): 
        v = list(self.vet0_l)
        x = list(self.vet0_l)
        a = list(self.vet0_l)
        for i in range(self.nb):
           xk = list(x0)
           v1 = v0[i]
           a1 = f(i, m, x0)
           v2 = v0[i] + a1 * dt / 2.
           xk[i] = x0[i] + v1 * dt / 2.
           a2 = f(i, m, xk)
           v3 = v0[i] + a2 * dt / 2.
           xk[i] = x0[i] + v2 * dt / 2.
           a3 = f(i, m, xk)
           v4 = v0[i] + dt * a3
           xk[i] = x0[i] + v3 * dt
           a4 = f(i, m, xk)
           x[i] = x0[i] + dt * (v1 + 2.*v2 + 2.*v3 + v4) / 6.0
           v[i] = v0[i] + dt * (a1 + 2.*a2 + 2.*a3 + a4) / 6.0
           a[i] = (a1 + 2.*a2 + 2.*a3 + a4) / 6.0
        return [x, v, a]
    
    # processes all steps
    def calculate(self, m, x0, v0, dt, T, method):
        self.nb = len(m)  # number of bodies
        t = x = v = a = []
        t.append(0)
        self.vet0_l = [array([.0, .0, .0]) for i in range(self.nb)]
        x = [list(self.vet0_l)]
        v = [list(self.vet0_l)]
        a = [list(self.vet0_l)]
        t[0] = 0
        x[0] = x0
        v[0] = v0
    
        for i in range(int(T / dt) - 1):
           t.append(0)
           x.append(list(self.vet0_l))
           v.append(list(self.vet0_l))
           a.append(list(self.vet0_l))
    
           if(method == 'rk4'):
              [x[i + 1], v[i + 1], a[i + 1]] = self.rk4(self.f, m, x[i], v[i], dt)
           else:
              [x[i + 1], v[i + 1], a[i + 1]] = self.eulerm(self.f, m, x[i], v[i], dt)
    
           # time
           t[i + 1] = t[i] + dt
    
        return [t, x, v, a]


# x[i][n][c]
# i = vector index of all positions at the same time
# n = identifies the body
# c = identifies the component [x=0,y=1,z=2]
# x.pop(0) remove older

