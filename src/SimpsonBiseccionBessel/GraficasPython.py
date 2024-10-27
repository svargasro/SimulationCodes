import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sci

Nsteps=100

xdata=np.zeros(Nsteps)
ydata=np.zeros(Nsteps)

f = open('datos.dat', 'r')
i=0
for line in f:
    line = line.strip()
    columns = line.split()
    xdata[i] = float(columns[0])
    ydata[i] = float(columns[1])
    i+=1
f.close()


plt.plot(xdata,ydata)
plt.show()
