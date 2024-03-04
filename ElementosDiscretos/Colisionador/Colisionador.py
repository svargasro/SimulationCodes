import numpy as np
import matplotlib.pyplot as plt

#Constantes del problema físico.

N=2
G=1.0

#Constantes del algoritmo de integración
xi=0.1786178958448091
Lambda=-0.2123418310626054
chi=-0.06626458266981849
Um2lambdau2=(1-2*Lambda)/2
Um2chiplusxi=1-2*(chi+xi)

class Cuerpo:
    def __init__(self,x0,y0,z0,Vx0,Vy0,Vz0,m0,R0):
        self.m=m0
        self.r=np.array([x0,y0,z0])
        self.V=np.array([Vx0,Vy0,Vz0])
        self.R = R0
    def BorreFuerza(self):
        self.F = np.array((0,0,0))
    def SumeFuerza(self,dF):
        self.F = self.F + dF

    def Mueva_r(self,dt, coeficiente):
        self.r = self.r + self.V*(coeficiente*dt)

    def Mueva_V(self, dt, coeficiente):
        self.V= self.V+ self.F*(coeficiente*dt/self.m)

class Colisionador:
    def CalculeFuerzaEntre(self,Planeta1, Planeta2):
        m1=Planeta1.m
        m2=Planeta2.m
        r21 = Planeta2.r-Planeta1.r
        r2= np.sum(r21**2)
        aux = G*m2*m1*pow(r2,-1.5)
        F1 = r21*aux
        Planeta1.SumeFuerza(F1)
        Planeta2.SumeFuerza(-1*F1)

    def CalculeTodasLasFuerzas(self,Planeta):
        #Borro las fuerzas de todos los planeta
        for i in range(0,N):
            Planeta[i].BorreFuerza()
            #Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
        for i in range(0,N):
            for j in range(0,i):
                self.CalculeFuerzaEntre(Planeta[i],Planeta[j])




r=11
m0=10
m1=1
M=m0+m1
mu=m0*m1/M
x0=-m1*r/M
x1=m0*r/M
omega=np.sqrt(G*M/(r*r*r))
T=2*np.pi/omega
V0=omega*x0
V1=omega*x1
dt=0.1
ttotal=T

#INICIO
#---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)

Planeta0 = Cuerpo(x0, 0, 0,  0, 0.5*V0,  0,m0,1.0)
Planeta1 = Cuerpo(x1, 0, 0,  0, 0.5*V1,  0,m1,0.5)
Planeta = np.array((Planeta0, Planeta1))
Newton = Colisionador()


tArray = np.arange(0,ttotal + dt, dt)

infoPosicionPlaneta0x = np.zeros(tArray.shape)
infoPosicionPlaneta0y = np.zeros(tArray.shape)
infoPosicionPlaneta1x = np.zeros(tArray.shape)
infoPosicionPlaneta1y = np.zeros(tArray.shape)


#CORRO

for i,t in enumerate(tArray):

    infoPosicionPlaneta0x[i] = Planeta[0].r[0]
    infoPosicionPlaneta0y[i] = Planeta[0].r[1]
    infoPosicionPlaneta1x[i] = Planeta[1].r[0]
    infoPosicionPlaneta1y[i] = Planeta[1].r[1]

    for i in range(0,N):
        Planeta[i].Mueva_r(dt,xi)

    Newton.CalculeTodasLasFuerzas(Planeta)

    for i in range(0,N):
        Planeta[i].Mueva_V(dt,Um2lambdau2)

    for i in range(0,N):
        Planeta[i].Mueva_r(dt,chi)

    Newton.CalculeTodasLasFuerzas(Planeta)

    for i in range(0,N):
        Planeta[i].Mueva_V(dt,Lambda)

    for i in range(0,N):
        Planeta[i].Mueva_r(dt,Um2chiplusxi)

    Newton.CalculeTodasLasFuerzas(Planeta)

    for i in range(0,N):
        Planeta[i].Mueva_V(dt,Lambda)

    for i in range(0,N):
        Planeta[i].Mueva_r(dt,chi)

    Newton.CalculeTodasLasFuerzas(Planeta)

    for i in range(0,N):
        Planeta[i].Mueva_V(dt,Um2lambdau2)

    for i in range(0,N):
        Planeta[i].Mueva_r(dt,xi);

plt.style.use('seaborn-v0_8')

fig, axes = plt.subplots(2,1,figsize=(6, 7))


axes[0].plot(infoPosicionPlaneta0x, infoPosicionPlaneta0y, '-', markersize=2, label=r'Planeta 1')

# Se ajustan demás detalles del gráfico.
axes[0].set_xlabel('x(t)', fontsize=12)
axes[0].set_ylabel('y(t)', fontsize=12)
# axes.legend(loc='upper left')
axes[0].grid(True, linestyle='--')
axes[0].set_title("Planeta 1.", fontsize=14)


axes[1].plot(infoPosicionPlaneta1x, infoPosicionPlaneta1y, '-', markersize=2, label=r'Planeta 2')


# Se ajustan demás detalles del gráfico.
axes[1].set_xlabel('x(t)', fontsize=12)
axes[1].set_ylabel('y(t)', fontsize=12)
# axes.legend(loc='upper left')
axes[1].grid(True, linestyle='--')
axes[1].set_title("Planeta 2.", fontsize=14)

plt.tight_layout()
plt.show()
fig.savefig('Colisionador.pdf')
