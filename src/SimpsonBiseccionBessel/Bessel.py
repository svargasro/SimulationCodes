def f(t,x,n):
  return np.cos(n*t-x*np.sin(t))



def IntegralPorSimpson(a,b,nsteps,x,n):
  nsteps*=2
  h=(b-a)/nsteps
  suma=0
  for i in range(nsteps+1):
    t=a+i*h
    if(i==0 or i==nsteps):
      suma+=f(t,x,n)
    elif(i%2==0):
      suma+=2*f(t,x,n)
    else:
      suma+=4*f(t,x,n)
  return suma*h/3


def Bessel(n,x):
  return (1/np.pi)*IntegralPorSimpson(0,np.pi,100,x,n)
