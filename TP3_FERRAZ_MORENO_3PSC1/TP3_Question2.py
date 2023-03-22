# -*- coding: utf-8 -*-
"""
TP 3 

FERRAZ Marc-Aurele
MORENO Antonin

3PSC1
"""
import numpy as np
from math import sin
from math import cos
from math import exp
from math import sqrt
import matplotlib.pyplot as plt

def VectEulerExp(f, a, b, ic, N):
    h = (b - a) / N #step size if h is constant      
    Lt = np.linspace(a, b, N) 
    Ly = np.empty((N, np.size(ic)), dtype = float) 
    Ly[0,:] = ic      
    for i in range(N-1):
          #if h isn't constant, we use h=t[i+1]-t[i]          
          Ly[i+1,:] = Ly[i,:] + h*f(Lt[i],Ly[i,:])
          # yn+1 = yn  +h f(tn,yn)      
    return (Lt, Ly)

def ode_RK2(f, a, b, ic, N):
    h = (b - a) / N #step size if h is constant
    Lt = np.linspace(a, b, N)
    Ly = np.empty((N, np.size(ic)),dtype = float)
    Ly[0,:] = ic
    for i in range(N-1):
    #if h isn't constant, we use h = t[i+1]-t[i]
        k1 = h*f(Lt[i], Ly[i,:])
        k2 = h*f(Lt[i] + h/2, Ly[i,:]+k1/2)
        Ly[i+1,:] = Ly[i,:] + k2
    return (Lt, Ly)


def ode_RK4(f, a, b,ic, N):
    h = (b - a) / N #step size if h is constant
    Lt = np.linspace(a, b, N)
    Ly = np.empty((N, np.size(ic)),dtype = float)
    Ly[0,:] = ic
    for i in range(N-1):
    #if h isn't constant, we use h=t[i+1]-t[i]
        k1 = h*f(Lt[i], Ly[i,:])
        y1 = Ly[i,:] + 1/2*k1
        k2 = h* f(Lt[i]+h/2, y1)
        y2 = Ly[i,:] + 1/2*k2
        k3 = h* f(Lt[i]+h/2,y2)
        y3 = Ly[i,:] + k3
        k4 = h* f(Lt[i]+h, y3)
        k = (k1+2*k2+2*k3+k4)/6
        Ly[i+1,:] = Ly[i,:] + k
    return (Lt, Ly)

def f2(t,Y):
    [y,dy] = Y
    return np.array([dy, sin(y) + sin(t)])


a1 = 0
b1 = 10
ic2 = np.array([0,-1]) 
N = 50
plt.figure(1)
#Test of ode_EulerExp == Runge-Kutta 1
T,Y= VectEulerExp(f2, a1, b1, ic2, N)
Y = Y[:,0]
plt.plot(T,Y,label = 'Explicit Euler')
T,Y= ode_RK2(f2, a1, b1, ic2, N)
Y = Y[:,0]
plt.plot(T,Y,label = 'RK2')
T,Y= ode_RK4(f2, a1, b1, ic2, N)
Y = Y[:,0]
plt.plot(T,Y,label = 'RK4')


#Y  = [f2_exact(t) for t in T]
#plt.plot(T,Y, label = 'Theoretical solution')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title("Comparaison des méthode de résolution différentielle numérique n = 50")
plt.legend()
plt.grid()
plt.show()