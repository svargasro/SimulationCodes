//CA de Difusion continuo 1D en C++
#include  <iostream>
#include <fstream>
#include  <cmath>
#include "Random64.h"
using namespace std;



/////////----- CONSTANTES GLOBALES
const int Lx=128;
const int Ly=128;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5; //Tiempo de relajación.
const double Utau=1.0/tau;
const double UmUtau=1-Utau;



class LatticeBoltzmann{
private:
double w[Q];      //Weights
  int Vx[Q],Vy[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions

  public:
    LatticeBoltzmann(void);
    ~LatticeBoltzmann(void);
    int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
    double rho(int ix,int iy,bool UseNew);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}


double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
//------------------- FUNCIONES GLOBALES -------

int main(void){

  
  return 0;
}
