import numpy as np
import matplotlib.pyplot as plt

ErrMax=1e-7

def f (x):
    return np.sin(x)/x

def CerosPorBiseccion(a,b):
    fa=f(a)

    while (b-a)>ErrMax:
        m = (a+b)/2.0
        fm = f(m)

        if (fm*fa<0):
            b=m
        else:
            a=m
            fa=fm

    return (a+b)/2


print("El cero es x= ",CerosPorBiseccion(2,4))
