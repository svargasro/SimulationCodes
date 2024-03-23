#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

const double Gamma = 16.0;
const double Deltat =0.001;
const double alpha = 1- exp(-Gamma*Deltat);
const double alphaAux = alpha*(2-alpha);
const int N = 1000;

//Deaclaración de la clase
class Cuerpo;

//Deaclaración de la interfase
class Cuerpo{
private:
double  x,Vx,Fx; double m,R;
public:
  void Inicie(double x0, double Vx0, double m0,double R0);
  void CalculeFuerza(void);
  void Arranque(double dt);
  void Muevase(double dt, double kT, Crandom & ran64);
  double Getx(void){return x;}; // Inline
  double GetVx(void){return Vx;}; // Inline
};
//Implementación de las funciones
void Cuerpo::Inicie(double x0, double Vx0, double m0,double R0){
  x=x0; Vx= Vx0; m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
Fx=0;
}

//Algoritmo de LeapFrog
void Cuerpo::Arranque(double dt){
  Vx=-(dt/(2*m))*Fx;
}

void Cuerpo::Muevase(double dt, double kT, Crandom & ran64){
  //Algoritmo Browniano de Leap Frog.
  double Vprime = Vx +Fx*dt;
  double DeltaV = -alpha*Vprime+sqrt(alphaAux*kT/m)*ran64.gauss(0,1);
  x += (Vprime + DeltaV/2.0)*dt;
  Vx = Vprime + DeltaV;

}

//----------- Funciones Globales -----------

double Sigma2(Cuerpo * Polen){
  double suma, xprom, sigma2;
  int i;
  for(suma=0; i<N; i++){
    suma+= Polen[i].Getx();
  }
  xprom = suma/N;

    for(suma=0; i<N; i++){
    suma+= pow (Polen[i].Getx() - xprom,2.0);
  }
    sigma2 = suma/(N-1);
  return sigma2;

}


int main(){


  Cuerpo Polen[N];

  Crandom ran64(1);


  double r0=1;
  double m0=1;
  double kT=4.0;
  double t,ttotal=1;
  int Ncuadros=500;
  double tdibujo,tcuadro=ttotal/Ncuadros;



  for(int i=0; i<N; i++) Polen[i].Arranque(Deltat);

  //----------(x0,Vx0,m0 ,R0)
  for(int i=0; i<N; i++) Polen[i].Inicie(0,0,m0,r0);
  
  
    //   cout<<t<<" "<<Polen.Getx()<<" "<<Polen.GetVx()<<endl;
   for(t=tdibujo=0;t<ttotal;t+=Deltat,tdibujo+=Deltat){
     if(tdibujo>tcuadro){
    //   tdibujo=0;
     cout<<t<<Sigma2(Polen)<<endl;
   }
    
    for(int i=0; i<N; i++) Polen[i].CalculeFuerza();
    for(int i=0; i<N; i++) Polen[i].Muevase(Deltat, kT, ran64);

  }
//  for(int i=0; i<N; i++) cout<<Polen[i].Getx()<<endl;
  return 0;
}
