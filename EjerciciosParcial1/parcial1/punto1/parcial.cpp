#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes del problema físico
const int N=1;
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
  double Getz(void){return r.z();}; // Inline
  friend class Colisionador;
};
class Colisionador{
private:
public:
  void CalculeTodasLasFuerzas(Cuerpo &balon);
  void CalculeFuerzaEntre(Cuerpo & balon1,Cuerpo & balon2);
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
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.z()<<"+"<<R<<"*sin(t)";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo &balon){
  int i,j;
  //Borro las fuerzas de todos los planetas
  balon.BorreFuerza();
  double g=9.8;
  vector3D Fg;
  Fg.load(0,0,-balon.m*g);
  balon.SumeFuerza(Fg);

  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  // for(i=0;i<N;i++)
  //   for(j=0;j<i;j++)
  //     CalculeFuerzaEntre(balon[i],balon[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & balon1,Cuerpo & balon2){
  double m1=balon1.m, m2=balon2.m;
  vector3D r21=balon2.r-balon1.r; double r2=r21.norm2();
  double aux=G*m2*m1*pow(r2,-1.5);
  vector3D F1=r21*aux;
  balon1.SumeFuerza(F1);  balon2.SumeFuerza(F1*(-1));
}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'balon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1:30]"<<endl;
  cout<<"set yrange[0:5]"<<endl;
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

  double m0 = 0.43;
  double x0=0, y0=0, z0=0;
  double theta0 = 20.0*M_PI/180.0; //gamma
  double v0= 20;
  double vx0 = v0*cos(theta0), vy0=0, vz0=v0*sin(theta0);
  double r0= 0.22;
  double t,dt=1e-4,ttotal= 2*vz0/9.8;
  int Ncuadros=100; double tdibujo,tcuadro = ttotal/Ncuadros;
  Cuerpo balon;
  Colisionador Newton;
  int i;

//  double alcanceXTeorico = vx0*ttotal;

  InicieAnimacion();
  
  //INICIO
  //---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)
  balon.Inicie(x0, y0, z0,  vx0, vy0, vz0,m0,r0);

  clog<<dt<<"\t"<<dt<<endl; //Se pasa dt al archivo.
  //CORRO
  for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){

      InicieCuadro();
      balon.Dibujese();
      TermineCuadro();

      tdibujo=0;
    }
    clog<<balon.Getx()<<" "<<balon.Getz()<<endl;


  balon.Mueva_r(dt,xi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,Um2lambdau2);
  balon.Mueva_r(dt,chi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,lambda);
  balon.Mueva_r(dt,Um2chiplusxi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,lambda);
  balon.Mueva_r(dt,chi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,Um2lambdau2);
  balon.Mueva_r(dt,xi);

  }



return 0;
}
