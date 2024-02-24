import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

g=9.8
Deltat=0.001
Nsteps= 50
G=M=1.0


class Cuerpo:
    #Constructor
    def __init__(self, x0,y0,z0,Vx0,Vy0,Vz0,m0):
        self.r=np.array((x0,y0,z0))
        self.m=m0
        self.V=np.array((Vx0,Vy0,Vz0))

    def CalculeFuerza(self):

        self.F = (-G*M*self.m/np.linalg.norm(self.r)**3)*self.r

    def Muevase(self,dt):
       self.r=self.r+dt*self.V
       self.V=self.V+(dt/self.m)*self.F

r0=10
omega = np.sqrt(G*M/(r0**3))
T=2*np.pi/omega
V0 = omega*r0
xdata=np.zeros(Nsteps)
ydata=np.zeros(Nsteps)
Balon = Cuerpo(x0=r0, y0=0, z0=0, Vx0=0, Vy0=0.5*V0, Vz0=0, m0=0.453)

for i in range(Nsteps):
    t=i*Deltat
    xdata[i] = Balon.r[0]
    ydata[i] = Balon.r[1]
    #print(xdata[i], " ", ydata[i])
    Balon.CalculeFuerza()
    Balon.Muevase(Deltat)

#print(xdata, " ", ydata)

#plt.plot(xdata,ydata)
#plt.show()
