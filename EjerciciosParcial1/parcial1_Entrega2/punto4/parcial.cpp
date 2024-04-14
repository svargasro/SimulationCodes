#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "vector.h"
using namespace std;


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
  vector3D r,V,F,omega; double m,R,rho,theta,tau,I;
public:
  void Inicie(double x0,double y0,double z0, double Vx0,double Vy0,double Vz0,double m0,double R0, double rho0, double omega0, double theta0);
  void BorreFuerza(void){F.load(0,0,0); tau=0;};// Inline
  void SumeFuerza(vector3D dF, double dtau){F+=dF; tau+=dtau;};// Inline

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
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0,double Vz0,double m0,double R0, double rho0, double omega0, double theta0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0; rho=rho0;
  omega.load(0,0,omega0);
  theta = theta0;
  I=2.0/(5.0*m*R*R);
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt); theta+= omega.z()*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
  double omegaAnt = omega.z();
  omega.load(0,0, omegaAnt + tau*(coeficiente*dt/I));
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.z()<<"+"<<R<<"*sin(t)";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo &balon){
  //Borro las fuerzas de todos los planetas
  balon.BorreFuerza();


  double g=9.8;
  vector3D Fg;
  Fg.load(0,0,-balon.m*g);
  balon.SumeFuerza(Fg,0); //Fuerza de gravedad.

  vector3D vBalon = balon.V;
  double vNorm = vBalon.norm();
  double rhoB = balon.rho;
  double radioB = balon.R;
  double areaB = M_PI*pow(radioB,2);
  double cA = 0.5;
  vector3D Farrastre = vBalon*(-0.5)*cA*rhoB*areaB*vNorm;
  balon.SumeFuerza(Farrastre,0);

  double cM = 1.0;
  vector3D omegaBalon = balon.omega;
  double angleZita = angle(omegaBalon,vBalon); //De la librería vector.h
  vector3D omegaCrossV = omegaBalon^vBalon;
  vector3D FMagnus = -(0.5)*cM*rhoB*areaB*radioB*(1.0/sin(angleZita))*omegaCrossV;

  double betta = 2e-4;

  vector3D torqueMagnus = -betta*vNorm*omegaBalon;
  double torquez = torqueMagnus.z();

  balon.SumeFuerza(FMagnus,torquez);

}

//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
   // cout<<"set terminal gif animate"<<endl;
   // cout<<"set output 'balon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1:14]"<<endl;
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

 const std::vector<double> frecuencies {0.1,0.2,0.4,0.8,1.6,2.0,3.2};

 ofstream fout;
 fout.open("TablaydesvVSF_ConTorque.txt");

 fout<<"frecuencia"<<"\t"<<"xdesv"<<"\t"<<"ydesv"<<"\n";

  for (double frecuency : frecuencies){

  double m0 = 0.43;
  double x0=0, y0=0, z0=0;
  double theta0 = 20.0*M_PI/180.0; //gamma
  double v0= 20;
  double vx0 = v0*cos(theta0), vy0=0, vz0=v0*sin(theta0);
  double r0= 0.22;
  double rho0 = 1.224;
  double frec = frecuency;
  double omega0z = frec*2*M_PI;
  double theta0Balon = 0;

  double t,dt=1e-4,ttotal= 2*vz0/9.8; //Tiempo de vuelo sin fricción del aire. (Se verifica nuevamente 1e-4)
  int Ncuadros=1000; double tdibujo,tcuadro = ttotal/Ncuadros;
  Cuerpo balon;
  Colisionador Newton;




//  InicieAnimacion();
  
  //INICIO
  //---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)
  balon.Inicie(x0, y0, z0,  vx0, vy0, vz0,m0,r0,rho0, omega0z,theta0Balon);

//  clog<<dt<<"\t"<<dt<<endl; //Se pasa dt al archivo.


  //CORRO

  for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){
/*
  if(tdibujo>tcuadro){

    InicieCuadro();
    balon.Dibujese();
    TermineCuadro();

    tdibujo=0;
    }
*/


    double y = balon.Gety();
    double x = balon.Getx();

    if(fabs(1 - x/9.0) < 1e-4) {fout<<frec<<" "<<x<<" "<<y<<endl; break;}




  balon.Mueva_r(dt,xi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,Um2lambdau2);
  balon.Mueva_r(dt,chi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,lambda);
  balon.Mueva_r(dt,Um2chiplusxi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,lambda);
  balon.Mueva_r(dt,chi);
  Newton.CalculeTodasLasFuerzas(balon); balon.Mueva_V(dt,Um2lambdau2);
  balon.Mueva_r(dt,xi);


  if((balon.Getz())< 1e-3 && (t>ttotal/2.0)) break;
  }



  }

  fout.close();

return 0;
}
