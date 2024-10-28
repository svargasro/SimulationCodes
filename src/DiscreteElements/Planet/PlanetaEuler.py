import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sci
GM=1.0


class Cuerpo:
  def __init__(self,x0,y0,z0,Vx0,Vy0,Vz0,m0):
    self.m=m0
    self.r=np.array([x0,y0,z0])
    self.V=np.array([Vx0,Vy0,Vz0])
  def CalculeFuerza(self):
    aux=-GM*self.m/np.linalg.norm(self.r)**3
    self.F = aux*self.r
  def Muevase(self,dt):
    self.r=self.r+dt*self.V
    self.V=self.V+(dt/self.m)*self.F

Deltat=0.001
Nsteps=200000


r0=10
omega=np.sqrt(GM/(r0**3))
T=2*np.pi/omega
V0=omega*r0
xdata=np.zeros(Nsteps)
ydata=np.zeros(Nsteps)
Balon = Cuerpo(x0=r0,y0=0,z0=0,Vx0=0,Vy0=0.5*V0,Vz0=0,m0=0.453) # Creo el ejemplar Balon de la clase Cuerpo
for i in range (Nsteps):
  t=i*Deltat
  xdata[i] = Balon.r[0]
  ydata[i] = Balon.r[1]
  Balon.CalculeFuerza()
  Balon.Muevase(Deltat)

plt.plot(xdata,ydata)
plt.show()
