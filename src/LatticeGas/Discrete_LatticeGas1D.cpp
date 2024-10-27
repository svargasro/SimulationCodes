//CA de Difusion 1D en C++
#include  <iostream>
#include  <cmath>
#include "Random64.h"
using namespace std;

const int Lx=1024;
const double p=0.5;

const int Q=2;

class LatticeGas{
private:
  int V[Q]; //V[i] i=0 (derecha) i=1 (izquierda)
  int n[Lx][Q],nnew[Lx][Q]; // n[ix][i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double mu,double sigma,Crandom & ran64);
  double rho(int ix,bool UseNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  void Show(void);
  double GetSigma2(void);
};
LatticeGas::LatticeGas(void){
  //definir los vectores velocidad
  V[0]=1;  V[1]=-1;
}
void LatticeGas::Borrese(void){
  for(int ix=0;ix<Lx;ix++)
    for(int i=0;i<Q;i++)
      n[ix][i]=0;
}
void LatticeGas::Inicie(int N,double mu,double sigma,Crandom & ran64){
  int ix,i;
  while(N>0){
    //Escojo un sitio al azar usando una distribucion gaussiana;
    ix=(int) ran64.gauss(mu,sigma); if(ix<0) ix=0; if(ix>=Lx) ix=Lx-1;
    //Escojo al azar entre las dos direcciones;
    i=(int) 2*ran64.r();
    //Si ese sitio esta vacio, lo lleno y decremento n en 1;
    if(n[ix][i]==0) {n[ix][i]++; N--;}
  }
}
double LatticeGas::rho(int ix,bool UseNew){
  if(UseNew)
    return nnew[ix][0]+nnew[ix][1];
  else
    return n[ix][0]+n[ix][1];
}

void LatticeGas::Colisione(Crandom & ran64){
  int ix,i;
  for(ix=0;ix<Lx;ix++) //para cada celda
    if(ran64.r()<p) //Genero un número al azar, y si es menor que p
      for(i=0;i<Q;i++) nnew[ix][i]=n[ix][i]; //Lo dejo igual
    else
      for(i=0;i<Q;i++) nnew[ix][i]=n[ix][(i+1)%2];  //Intercambio los contenidos
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<Q;i++) // y en cada dirección
      n[(ix+Lx+V[i])%Lx][i]=nnew[ix][i]; //Condición de frontera periódica.
}
void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){
    for(int ix=0;ix<Lx;ix++)
      cout<<n[ix][i]<<" ";
    cout<<endl;
  }
  cout<<endl;
}
double LatticeGas::GetSigma2(void){
  int ix;
  //Calcular cuántas bolitas hay
  double N=0;
  for(ix=0;ix<Lx;ix++)
    N+=rho(ix,false);
  //Calcular la posición promedio  
  double xprom=0;
  for(ix=0;ix<Lx;ix++)
    xprom+=ix*rho(ix,false);
  xprom/=N;
  //Calcular la varianza promedio
  double Sigma2=0;
  for(ix=0;ix<Lx;ix++)
    Sigma2+=pow(ix-xprom,2.0)*rho(ix,false);
  Sigma2/=(N-1);

  return Sigma2;
}

int main(void){

  LatticeGas Difusion;
  Crandom ran64(1);
  //int N= (int)Lx*0.1;
  int N = 128;
  double mu=Lx/2, sigma=Lx/8;
  int t,tmax=400;

  Difusion.Borrese();
  Difusion.Inicie(N,mu,sigma,ran64);

  for(t=0;t<tmax;t++){
    cout<<t<<" "<<Difusion.GetSigma2()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  
  
  return 0;
}
