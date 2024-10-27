import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

g=9.8
Deltat=0.01
Nsteps=200



class Cuerpo:
    #Constructor
    def __init__(self, x0,y0,z0,Vx0,Vy0,Vz0,m0):
        self.r=np.array((x0,y0,z0))
        self.m=m0
        self.V=np.array((Vx0,Vy0,Vz0))

    def CalculeFuerza(self):
        self.F=np.array((0,-self.m*g,0))

    def Muevase(self,dt):
       self.r=self.r+dt*self.V
       self.V=self.V+(dt/self.m)*self.F

xdata=np.zeros(Nsteps)
ydata=np.zeros(Nsteps)
Balon = Cuerpo(x0=0, y0=0, z0=0, Vx0=16, Vy0=9, Vz0=0, m0=0.453)

for i in range(Nsteps):
    t=i*Deltat
    xdata[i] = Balon.r[0]
    ydata[i] = Balon.r[1]
    Balon.CalculeFuerza()
    Balon.Muevase(Deltat)

#plt.plot(xdata, ydata)
#plt.show()

# Animación

# fig = plt.figure()
# ax = plt.axes(xlim=(0,40), ylim=(-8,4.5))
# line, = ax.plot([], [],'o')

# def init():
#     line.set_data([], [])
#     return line,
# def animate(i):
#     x=xdata[i]
#     y=ydata[i]
#     line.set_data(x, y)
#     return line,


# anim = animation.FuncAnimation(fig, animate, init_func=init,frames=Nsteps, interval=20, blit=False)
# anim.save('BalonAnimado.mp4', fps=30)
