import numpy as np
from math import *

import matplotlib.pyplot as plt


def method_S(f, a, b, n):
    h = (b-a)/n
    x = np.linspace(a, b-h, n)
    Sx = sum(f(x[1:]))
    Sm = sum(f(x+h/2))
    A = f(a) + f(b) + 2*Sx + 4*Sm
    return h*A/6


def method_RD(f, a, b, n):
    h = (b-a)/n
    x = np.linspace(a+h,b,n)
    S = sum(f(x))
    return h*S

def f(t):
    return t**(3/2)*np.exp(-t/2)

def F(x):
    return (1/(3*sqrt(2*pi))*method_RD(f, 0, x, n))

def F1(x):
    return (1/(3*sqrt(2*pi))*method_S(f, 0, x, n))







liste_y = []
liste_x = []
liste_y1 = []
liste_x1 = []

for x in np.arange(0,20.1,0.1):
    n=50
    liste_y.append(F(x))
    liste_x.append(x)
    liste_y1.append(F1(x))
    liste_x1.append(x)


print(liste_y[-1])
print(liste_y1[-1])
    
plt.plot(liste_x, liste_y)
plt.plot(liste_x1, liste_y1)
plt.xlabel("x")
plt.ylabel("F(x)")
plt.legend(["Méthode RD", "Méthode Simpson"], loc ="lower right")
plt.title('F en fonction de x, avec différente méthode de calcul intégral avec n = 50')
plt.show()

