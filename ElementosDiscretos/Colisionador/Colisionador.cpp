#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes del problema físico
const int N=2;
const double G=1.0;

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

//--------------- Declarar las clases-----------
class Cuerpo;
class Colisionador;

//--------- Declarar las interfases de las clases---------
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,
              double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline

  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  friend class Colisionador;
};
class Colisionador{
private:
public:
  void CalculeTodasLasFuerzas(Cuerpo * Planetas);
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Planeta){
  int i,j;
  //Borro las fuerzas de todos los planetas
  for(i=0;i<N;i++)
    Planeta[i]. BorreFuerza();
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=0;j<i;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  double m1=Planeta1.m, m2=Planeta2.m;
  vector3D r21=Planeta2.r-Planeta1.r; double r2=r21.norm2();
  double aux=G*m2*m1*pow(r2,-1.5);
  vector3D F1=r21*aux;
  Planeta1.SumeFuerza(F1);  Planeta2.SumeFuerza(F1*(-1));
}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
   // cout<<"set terminal gif animate"<<endl;
   // cout<<"set output 'Colisionador.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-11:11]"<<endl;
  cout<<"set yrange[-11:11]"<<endl;
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

int main(){
  double r=11,m0=10,m1=1;
  double M=m0+m1, mu=m0*m1/M;
  double x0=-m1*r/M,x1=m0*r/M;
  double omega=sqrt(G*M/(r*r*r)); double T=2*M_PI/omega;
  double V0=omega*x0, V1=omega*x1;
  double t,dt=0.01,ttotal=T;

  int Ncuadros=30000; double tdibujo,tcuadro = ttotal/Ncuadros;
  Cuerpo Planeta[N];
  Colisionador Newton;
  int i;

 InicieAnimacion();
  
  //INICIO
  //---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)
  Planeta[0].Inicie(x0, 0, 0,  0, 0.5*V0,  0,m0,1.0);
  Planeta[1].Inicie(x1, 0, 0,  0, 0.5*V1,  0,m1,0.5);
  //CORRO
  cout<<tcuadro<<endl;
  for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){

      InicieCuadro();
      for(i=0;i<N;i++) Planeta[i].Dibujese();
      TermineCuadro();

      tdibujo=0;
    }
//  cout<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<endl;

  for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,xi);
  Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Um2lambdau2);
  for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,chi);
  Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,lambda);
  for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Um2chiplusxi);
  Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++)Planeta[i].Mueva_V(dt,lambda);
  for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,chi);
  Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++)Planeta[i].Mueva_V(dt,Um2lambdau2);
  for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,xi);
    
}
return 0;
}
