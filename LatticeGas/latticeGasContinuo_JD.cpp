//CA de Difusion continuo 1D en C++
#include  <iostream>
#include  <cmath>
#include "Random64.h"
using namespace std;

const int Lx=128;
const double p=0.5;

const int Q=2;

class LatticeGas{
private:
  int V[Q]; //V[i] i=0 (derecha) i=1 (izquierda)
  double f[Lx][Q],fnew[Lx][Q]; // f[ix][i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(double mu,double sigma);
  double rho(int ix,bool UseNew);
  void Colisione(void);
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
      f[ix][i]=0;
}
void LatticeGas::Inicie(double mu,double sigma){
  for(int ix=0;ix<Lx;ix++)
    for(int i=0;i<Q;i++)
      f[ix][i]=fnew[ix][i]=0.5/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2.0));
}
double LatticeGas::rho(int ix,bool UseNew){
  if(UseNew)
    return fnew[ix][0]+fnew[ix][1];
  else
    return f[ix][0]+f[ix][1];
}

void LatticeGas::Colisione(void){
  int ix,i,j;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<Q;i++)
      {j=(i+1)%2; fnew[ix][i]=p*f[ix][i]+(1-p)*f[ix][j];}
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<Q;i++) // y en cada dirección
      f[(ix+Lx+V[i])%Lx][i]=fnew[ix][i];
}
void LatticeGas::Show(void){
  for(int ix=0;ix<Lx;ix++)
    cout<<ix<<" "<<rho(ix,true)<<endl;
  cout<<endl;
}

double LatticeGas::GetSigma2(void){
  int ix;
  //Calcular cuántas bolitas hay
  double N=0;
  for(ix=0;ix<Lx;ix++)
      N+=rho(ix,true);
  //Calcular la posición promedio  
  double xprom=0;
  for(ix=0;ix<Lx;ix++)
    xprom+=ix*rho(ix,true);
  xprom/=N;
  //Calcular la varianza promedio
  double Sigma2=0;
  for(ix=0;ix<Lx;ix++)
    Sigma2+=pow(ix-xprom,2.0)*rho(ix,true);
  Sigma2/=N;

  return Sigma2;
}

//------------------- FUNCIONES GLOBALES -------


int main(void){
  LatticeGas Difusion;
  int N=1; double mu=Lx/2, sigma=Lx/64;
  int t,tmax=400;

  Difusion.Borrese();
  Difusion.Inicie(mu,sigma);

  for(t=0;t<tmax;t++){
    cout<<t<<" "<<Difusion.GetSigma2()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }

  //  Difusion.Show();
  
  
  return 0;
}
