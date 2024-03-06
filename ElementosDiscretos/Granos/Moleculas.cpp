#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//Constantes del problema físico

const double G=1.0;

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);
const double kHertz=1.0e4;
const int Nx=5;
const int Ny=5;
const int N=Nx*Ny;
const double Lx=60;
const double Ly=60;

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
  void CalculeTodasLasFuerzas(Cuerpo * Moleculas);
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
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
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
  //Borro las fuerzas de todos los planetas
  for(i=0;i<N;i++)
    Molecula[i]. BorreFuerza();
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=0;j<i;j++)
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  vector3D r21=Molecula2.r-Molecula1.r; double d=r21.norm();
  //Calculo si hay colisión.
  double R1 = Molecula1.R;
  double R2 = Molecula2.R;
  //Si d=R1+R2 se besan, si d> R1+R2 no se chocan, si d<R1+R2 entonces se están chocando.
  //La interperpenetración es S.
  double s= (R1+R2)-d;
  if (s>0) {//Si hay colisión, calculo el vector normal.
  vector3D n = r21*(1.0/d);
  //Calculo Fuerza
  vector3D F2=n*(kHertz*pow(s,1.5));
  //La sumo a los granos.
  Molecula1.SumeFuerza(F2);  Molecula2.SumeFuerza(F2 *(-1));
  }


}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-11:"<<Ly+10<<"]"<<endl;
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

  Cuerpo Molecula[N];
  Crandom ran64(1);
  Colisionador Newton;
  int ix,jy;


  //Parámetros caja:

  //Parámetros Simulación:
  double m0=1.0;
  double R0=2.0;
  double  kT=10;
  double V0 = sqrt(kT/m0);

  //Variables auxiliares.
  double dx= Lx/(Nx+1);
  double dy = Ly/(Ny+1);
  double x0,y0,Vx0,Vy0;

  int Ncuadros= 1000;
  double t,dt=1e-3,ttotal = 10*(Lx/V0), tdibujo,tcuadro=ttotal/Ncuadros ;

  InicieAnimacion();
  
  //INICIO
  for(ix=0;ix<Nx;ix++){
for (jy=0;jy<Ny ;jy++ ) {
  double  theta = 2*M_PI*ran64.r();
  x0=(ix+1)*dx; y0=(jy+1)*dy;Vx0=V0*cos(theta); Vy0=V0*sin(theta);
  //---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)

  Molecula[jy*Nx+ix].Inicie(x0, y0, 0,  Vx0, Vy0,  0,m0,R0);


 }
  }



  //CORRO
  for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){

      InicieCuadro();
      for(int i=0;i<N;i++) Molecula[i].Dibujese();
      TermineCuadro();

      tdibujo=0;
    }
  // cout<<Molecula[0].Getx()<<" "<<Molecula[0].Gety()<<endl;

  // for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,xi);
  // Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Um2lambdau2);
  // for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,chi);
  // Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,lambda);
  // for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Um2chiplusxi);
  // Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
  // for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,chi);
  // Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,Um2lambdau2);
  // for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,xi);
    
}
return 0;
}
