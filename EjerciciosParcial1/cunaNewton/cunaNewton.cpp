#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"

using namespace std;

//Constantes del problema físico
double Lx=9, Ly=12;
const int Nx=3, Ny=1;
const int N=Nx*Ny;
const double G=1.0;
const double KHertz=20e10;
const double g=9.8;
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

  public:
    vector3D r,V,F; double m,R,xp,yp;
    void Inicie(double x0,double y0,double z0,
                double Vx0,double Vy0,double Vz0,double m0,double R0, double Xp, double Yp);
    void BorreFuerza(void){F.load(0,0,0);};// Inline
    void SumeFuerza(vector3D dF){F+=dF;};// Inline
    void Mueva_r(double dt,double coeficiente);
    void Mueva_V(double dt,double coeficiente);
    void Dibujese(void);
    double Getx(void){return r.x();}; // Inline
    double Gety(void){return r.y();}; // Inline
    double GetVx(void){return V.x();}                                  //
    double GetVy(void){return V.y();}                                  //

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
                    double Vx0,double Vy0,double Vz0,double m0,double R0, double Xp, double Yp){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
  xp = Xp; yp = Yp;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  double x=r.x();
  double y = r.y();
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
  cout<<" , (1- (t/7))*"<<x<<"+ (t/7)*"<<xp<<", (1- (t/7))*"<<y<<"+ (t/7)*"<<yp;

}



//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
  //Borro las fuerzas de todos los planetas
  for(i=0;i<N;i++)
    Molecula[i].BorreFuerza();


  for (int i=0;i<N ;i++ ) {
    vector3D Fp;
    double x = Molecula[i].r.x();
    double y = Molecula[i].r.y();
    double v2 = Molecula[i].V.norm2();
    double m = Molecula[i].m;
    double xp = Molecula[i].xp;
    double l = Molecula[i].yp;
    double sinMol = (x-xp)/l;
    double cosMol = (l-y)/l;
    double T= m*v2/l + g*m*cosMol;
    Fp.load(-T*sinMol,-m*g+T*cosMol,0);
    Molecula[i].SumeFuerza(Fp);
  }

  // for(int i=0;i<N;i++){
  //   vector3D Fg;
  //   Fg.load(0,-Molecula[i].m*g,0);
  //   Molecula[i].SumeFuerza(Fg);
  //}

  //N=3
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      // Molecula[i].r;
     CalculeFuerzaEntre(Molecula[i],Molecula[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  //Determinar si hay colision
  vector3D r21=Molecula2.r-Molecula1.r; double d=r21.norm();
  double s=(Molecula1.R+Molecula2.R)-d;
  if(s>0){ //Si hay colisión
    //Calcular el vector normal
    vector3D n=r21*(1.0/d);
    //Calculo la fuerza
    vector3D F2=n*(KHertz*pow(s,1.5));
    //Las sumo a los granos
    Molecula2.SumeFuerza(F2);  Molecula1.SumeFuerza(F2*(-1));
  }
}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'cunaNewton.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-12e-2:15e-2]"<<endl;
  cout<<"set yrange[-4e-2:12e-2]"<<endl;
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
  Colisionador Newton;
  Crandom ran64(1);
  int i,ix,iy;
  //Parametros de la simulación
  double m0=100e-2; double R0=1.5e-2;
  double kT=10; 
  //Variables auxiliares para la condición inicial
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  double theta; double V0=sqrt(kT/m0);
  double x0,y0,Vx0,Vy0;

  double w0 = sqrt(g/12e-2);
  double T = 2*M_PI/w0;
  //Variables auxiliares para correr la simulacion
  int Ncuadros=200; double t,tdibujo,dt=1e-6,tmax=T,tcuadro=tmax/Ncuadros;
  //1000
  InicieAnimacion();
  
  //INICIO

  //Inicializo las paredes: 4 granos gigantes.


  //Inicializo los granos.

  // for(ix=0;ix<Nx;ix++)
  //   for(iy=0;iy<Ny;iy++){
  //     theta=2*M_PI*ran64.r();
  //     x0=(ix+1)*dx; y0=(iy+1)*dy; Vx0=V0*cos(theta); Vy0=V0*sin(theta);
  //     //----------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  //     Molecula[iy*Nx+ix].Inicie(x0,y0, 0,Vx0,Vy0,  0,m0,R0);
  //   }

  double l0 = 12e-2;

  double theta0 = 15.0*M_PI/180.0;
  double x00 = -l0*sin(theta0);
  double y00 = l0*(-cos(theta0)+1);
  double xp0 = 0;
  double yp0 = l0;





    //     //----------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0,xp,yp)
  Molecula[0].Inicie(x00,y00,0,0,0,0,m0,R0,xp0,yp0);

  for (int i=1;i<N;i++) {
    Molecula[i].Inicie((i*2*R0),0,0,0,0,0,m0,R0,(i*2*R0),l0);
  }



  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){
      
     InicieCuadro();
    for(i=0;i<N;i++) Molecula[i].Dibujese();
     // Molecula[2].Dibujese();
     TermineCuadro();
      
      tdibujo=0;
    }
  clog<<t<<" "<<pow(Molecula[0].Gety() - Molecula[0].yp,2)+pow(Molecula[0].Getx()-Molecula[0].xp,2)<<endl;

    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,xi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Um2chiplusxi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,xi);
  }

  return 0;
}
