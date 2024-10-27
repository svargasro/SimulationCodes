#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

const double g=9.8;

//Deaclaración de la clase
class Cuerpo;

//Deaclaración de la interfase
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
};
//Implementación de las funciones
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0){
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  F.load(0,-m*g,0);
}
void Cuerpo::Muevase(double dt){
  //Algoritmo de Euler
  r+=V*dt;  V+=F*(dt/m);
}

//----------- Funciones Globales -----------


int main(){
  double t,dt=0.001,ttotal=2.5;
//  int Ncuadros=200; double tdibujo,tcuadro=ttotal/Ncuadros;
  Cuerpo Balon;


  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0   ,R0)
  Balon.Inicie( 0, 0, 0, 16,  9,  0,0.453,1.0);

  //for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){

  for(t=0; t<ttotal; t+=dt){
    // if(tdibujo>tcuadro){

    //   tdibujo=0;
    // }

    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
    Balon.CalculeFuerza();
    Balon.Muevase(dt);

  }
  return 0;
}
