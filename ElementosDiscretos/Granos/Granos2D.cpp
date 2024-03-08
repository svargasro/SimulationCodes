#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//Constantes del problema físico
double Lx=60, Ly=60; 
const int Nx=5, Ny=5;
const int N=Nx*Ny;
const double g=9.8;
const double KHertz=1.0e4;
const double Gamma = 150;
const double Kcundall= 500, mu=0.4;

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
  vector3D r,V,F; double m,R, theta, omega, tau, I;
public:
  void Inicie(double x0,double y0,
	      double Vx0,double Vy0, double theta0, double omega0, double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0); tau=0;};// Inline
  void SumeFuerza(vector3D dF,double dtau){F+=dF; tau+=dtau;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  friend class Colisionador;
};
class Colisionador{
private:
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo * Planetas);
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,
	      double Vx0,double Vy0, double theta0, double omega0, double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
  theta=theta0;
  omega=omega0;
  I= (2.0/5.0)*m*R*R;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano){
  int i,j;
  //Borro las fuerzas de todos los planetas
  for(i=0;i<N;i++)
    Grano[i]. BorreFuerza();
  //sumo fuerza de gravedad.
  for (int i=0;i<N ;i++ ) {
  vector3D Fg;
  Fg.load(0,-Grano[i].m*g,0);
  Grano[i].SumeFuerza(Fg,0);
  }

  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeFuerzaEntre(Grano[i],Grano[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2){
  //Determinar si hay colision
  vector3D r21=Grano2.r-Grano1.r; double d=r21.norm();
  double R1= Grano1.R;
  double R2=Grano2.R;
  double s=(R1+R2)-d;
  if(s>0){ //Si hay colisión
    //Calcular el vector normal
    vector3D n=r21*(1.0/d);

    //Calculo la fuerza
    vector3D F2=n*(KHertz*pow(s,1.5));
    //Las sumo a los granos
    Grano2.SumeFuerza(F2,0);  Grano1.SumeFuerza(F2*(-1),0);
    //Calcule Fuerza Tangencial
    vector3D Rw; Rw.load(0,0,R2*Grano2.omega - R1*Grano1.omega);
    vector3D Vc = (Grano2.V - Grano1.V) - (Rw^n);
    vector3D Vcn = n*(Vc*n), Vct= Vc-Vcn;
    //Suma la fuerza de deformación plástica normal (Kuramoto-Kano)
    double m1 = Grano1.m; double m2=Grano2.m; double m12 = m1*m2/(m1+m2);
    F2 = (-Gamma*sqrt(s)*m12)*Vcn;
    Grano2.SumeFuerza(F2,0);  Grano1.SumeFuerza(F2*(-1),0);

  }
}

void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
}

//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha

}
void TermineCuadro(void){
    cout<<endl;
}

int main(){
  Cuerpo Grano[N];
  Colisionador Hertz;
  Crandom ran64(1);
  int i,ix,iy;
  //Parametros de la simulación
  double m0=1.0; double R0=2.0;
  double kT=10; 
  //Variables auxiliares para la condición inicial
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  double theta; double V0=sqrt(kT/m0);
  double x0,y0,Vx0,Vy0;
  //Variables auxiliares para correr la simulacion
  int Ncuadros=1000; double t,tdibujo,dt=1e-3,tmax=10*Lx/V0,tcuadro=tmax/Ncuadros;

  InicieAnimacion();
  
  //INICIO

  //Inicializo las paredes: 4 granos gigantes.
  //Variables auxiliares para las paredes.
  //
  double Rpared = 100*Lx;
  double Mpared = 100*m0;

  /// (x0,y0,Vx0,Vy0,m0,R0)
  Grano[N].Inicie(Lx/2,Ly+Rpared,0,0,0,0,Mpared,Rpared); //Pared arriba
  Grano[N+1].Inicie(Lx/2,-Rpared,0,0,0,0,Mpared,Rpared); //Pared abajo
  Grano[N+2].Inicie(Lx+Rpared,Ly/2,0,0,0,0,Mpared,Rpared); //Pared derecha
  Grano[N+3].Inicie(-Rpared,Ly/2,0,0,0,0,Mpared,Rpared); //Pared izquierda


  //Inicializo los granos.
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      theta=2*M_PI*ran64.r();
      x0=(ix+1)*dx; y0=(iy+1)*dy; Vx0=V0*cos(theta); Vy0=V0*sin(theta);
      //----------------(x0,y0,Vx0,Vy0,theta0,omega0, m0,R0)
      Grano[iy*Nx+ix].Inicie(x0,y0,Vx0,Vy0,  0,0,m0,R0);
    }
      
  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
    cout<<Grano[1].Getx()<<" "<<Grano[1].Gety()<<endl;

    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Um2chiplusxi);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);
  }
  return 0;
}
