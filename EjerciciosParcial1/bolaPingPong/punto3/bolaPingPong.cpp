#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"
#include "Random64.h"

using namespace std;


//Constantes del problema físico
double Lx=10, Ly=60;
const int Nx=1, Ny=1;
const int N=Nx*Ny, Ntot=N+1;
const double g=9.8, KHertz=1.0e4, Gamma=10;
const double ARaqueta = 1;
const double TRaqueta = 2.9;
const double omegaRaqueta = 2*M_PI/TRaqueta;

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
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0); tau=0;};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  void muevaseRaqueta(double t);
  friend class Colisionador;
};
class Colisionador{
private:
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];
public:
    void CalculeTodasLasFuerzas(Cuerpo * Grano,double dt);
    void Inicie(void);
    void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,double dt);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
  theta=theta0; omega=omega0; I=2.0/5.0*m*R*R;
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

void Cuerpo::muevaseRaqueta(double t){
  r.load(0, -100*Lx+ARaqueta*sin(omegaRaqueta*t),0);
}


//------- Funciones de la clase Colisionador --------
void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano,double dt){
  int i,j;
  //Borro las fuerzas de todos los granos
  for(i=0;i<N+4;i++)
    Grano[i].BorreFuerza();

  //sumo fuerza de gravedad
  vector3D Fg;
  for(i=0;i<N;i++){
    Fg.load(0,-Grano[i].m*g,0);
    Grano[i].SumeFuerza(Fg);
  }
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+1;j++)
//      Grano[i].r;
      CalculeFuerzaEntre(Grano[i],Grano[j],dt);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,double dt){
  //Grano 1 es pelota.
  //Determinar si hay colision
  vector3D r21=Grano2.r-Grano1.r; double d=r21.norm();
  double R1=Grano1.R,R2=Grano2.R;
  double s=(R1+R2)-d;

  if(s>0){ //Si hay colisión
    //Vectores unitarios
    vector3D n=r21*(1.0/d);

    //Fn (Hertz-Kuramoto-Kano)
    double m1=Grano1.m;
    double v1 = Grano1.V.y();
    double Fn=KHertz*pow(s,1.5)-Gamma*sqrt(s)*m1*v1;
    //Calcula y Cargue las fuerzas
    vector3D F1,F2;
    F2=n*Fn;
    F1=F2*(-1);
//    clog<<"F1: "<<F1.x()<<endl;

    if(F1.y()> 0) Grano1.SumeFuerza(F1);

  }
}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'pingPongCaos.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange["<<-Lx-2<<":"<<Lx+2<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}

void InicieCuadro(Cuerpo &granoPared){
    cout<<"plot 0,0 ";
    cout<<" , (1 - t/7.0)*"<<-Lx<<"+(t/7.0)*"<<Lx<<","<<granoPared.Gety()+100*Lx;        //pared de abajo

//    cout<<" , (1 - t/7.0)*"<<-Lx<<"+(t/7.0)*"<<Lx<<",0";        //pared de abajo
    // cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    // cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    // cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  Cuerpo Grano[Ntot];
  Colisionador Hertz;
  Crandom ran64(1);
  //Parametros de la simulación
  double m0=1.0; double R0=2.0;
  //Variables auxiliares para la condición inicial
  double x0,y0,Vx0,Vy0;
  //Variables auxiliares para las paredes
  double Rpared=100*Lx, Mpared=100*m0;
  //Variables auxiliares para correr la simulacion

  int Ncuadros=500; double t,tdibujo,dt=1e-2,tmax=200,tcuadro=tmax/Ncuadros;
//100000

 InicieAnimacion();
  
  //INICIO
  //Inicializar las paredes
  //------------------(  x0,    y0,Vx0,Vy0,theta0,omega0,    m0,    R0)
 Grano[N].Inicie(0,-Rpared,  0,  0,     0,     0, Mpared,Rpared); //Pared arriba



  x0 = 0;
  y0 = 30;
  Vx0 = 0;
  Vy0 = 0;

  Grano[0].Inicie(x0,y0,Vx0,Vy0, 0,0,m0,R0);

  int i;


  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){
      
      InicieCuadro(Grano[N]);
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }


    //clog<<t<<"\t"<<Grano[0].Gety()<<endl;

    clog<<t<<"\t"<<Grano[0].Gety()<<"\t"<<Grano[1].Gety()+100*Lx<<endl;

    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Um2chiplusxi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);
    Grano[N].muevaseRaqueta(t);
  }


  return 0;
}
