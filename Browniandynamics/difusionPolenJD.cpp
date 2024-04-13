#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

//Constantes Globales
const double gammaVal=16.0;
const double Deltat=0.001;
const double alpha=1-exp(-gammaVal*Deltat);
const double alphaAux=alpha*(2-alpha);
const int N=1000;

//Deaclaración de la clase
class Cuerpo;

//Deaclaración de la interfase
class Cuerpo{
private:
  double x,Vx,Fx; double m,R;
public:
  void Inicie(double x0,double Vx0,double m0,double R0);
  void CalculeFuerza(void);
  void Arranque(double dt);
  void Muevase(double dt,double kT,Crandom & ran64);
  void Dibujese(void);
  double Getx(void){return x;}; // Inline
  double GetVx(void){return Vx;}; // Inline
};
//Implementación de las funciones
void Cuerpo::Inicie(double x0,double Vx0,double m0,double R0){
  x=x0;  Vx=Vx0; m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  Fx=0;
}
void Cuerpo::Arranque(double dt){
  //Algoritmo LeapFrog
  Vx-=Fx*(0.5*dt/m);
}
void Cuerpo::Muevase(double dt,double kT,Crandom & ran64){
  //Algoritmo Browniano de Leap Frog (van Gunsteren)
  double Vprime=Vx+Fx*dt;
  double DeltaV=-alpha*Vprime+sqrt(alphaAux*kT/m)*ran64.gauss(0,1);
  x+=(Vprime+DeltaV/2)*dt;
  Vx=Vprime+DeltaV;
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<0<<"+"<<R<<"*sin(t)";
}

//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:10]"<<endl;
  cout<<"set yrange[-2:2]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}

double Sigma2(Cuerpo * Polen){
  double suma,xprom,sigma2;
  int i;
  //Calculo xprom
  for(suma=0,i=0;i<N;i++)
    suma+=Polen[i].Getx();
  xprom=suma/N;
  //Calculo sigma2
  for(suma=0,i=0;i<N;i++)
    suma+=pow(Polen[i].Getx()-xprom,2.0);
  sigma2=suma/(N-1);
  return sigma2;
}

int main(){
  Cuerpo Polen[N];
  Crandom ran64(1);
  double x0=0,Vx0=0,R0=1,m0=1;
  double kT=4.0;
  double t,ttotal=1;
  int i,Ncuadros=1000; double tdibujo,tcuadro=ttotal/Ncuadros;
  
  InicieAnimacion();
  
  //----------(x0,Vx0,m0,R0)
  for(i=0;i<N;i++) Polen[i].Inicie(x0,Vx0,m0,R0);

  for(t=tdibujo=0;t<ttotal;t+=Deltat,tdibujo+=Deltat){

    for(i=0;i<N;i++) Polen[i].CalculeFuerza();
    for(i=0;i<N;i++) Polen[i].Arranque(Deltat);
    if(tdibujo>tcuadro){
//      cout<<t<<" "<<Sigma2(Polen)<<endl;

      InicieCuadro();
      for(i=0;i<N;i++)Polen[i].Dibujese();
      TermineCuadro();

      tdibujo=0;
    }
    for(i=0;i<N;i++) Polen[i].CalculeFuerza();
    for(i=0;i<N;i++) Polen[i].Muevase(Deltat,kT,ran64);
  }
  //  for(i=0;i<N;i++) cout<<Polen[i].Getx()<<endl;
  
  return 0;
}
