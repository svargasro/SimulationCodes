#include  <iostream>
#include  <cmath>
#include "Random64.h"

using namespace std;

const int Lx=1024;
const double p=0.5;

const int Q=2;

class LatticeGas{
    private:
    int V[Q]; //V[i] i=0 (Derecha) i=1 (Izquierda)
    int n[Lx][Q], nnew[Lx][Q]; //n[ix][i]


    public:
        LatticeGas(void);
        void Borrese(void);
        void Inicie(int N,double mu, double sigma, Crandom &ran64);
        void Show(void);
        void Colisione(Crandom &ran64);
        void Adveccione(void);
        double rho(int ix,bool UseNew);
        double GetSigma2(void);
};

void LatticeGas::Borrese(void){
    for(int ix=0; ix<Lx; ix++){
        for(int i=0;i<Q; i++){
        n[ix][i]=0;
        }
    }
}


LatticeGas::LatticeGas(void){
    //definir los vectores velocidad.
    V[0] = 1;
    V[1] = -1;

}
void LatticeGas::Inicie(int N,double mu, double sigma, Crandom &ran64){

    int ix,i;

    while(N>0){
        //Se elige un sitio al azar usando una distribución gaussiana.
        ix= (int)ran64.gauss(mu,sigma); if(ix<0) ix=0; if(ix>=Lx) ix=Lx-1;

        //Escoger al azar entre las dos direcciones
        i = (int)2*ran64.r();
        //Si está vacío, lo lleno y decremento n en 1
        if(n[ix][i] == 0){n[ix][i]++; N--;}
    }
}

void LatticeGas::Show(void){

    for(int i=0;i<Q; i++){
    for(int ix=0; ix<Lx; ix++)
         cout<<n[ix][i]<<" ";

        cout<<endl;
    }
}


void LatticeGas::Colisione(Crandom &ran64){

    for(int ix=0;ix<Lx;ix++){
        if(ran64.r()<p)
            for(int i=0;i<Q;i++)
                nnew[ix][i] = n[ix][i];
        else
            for(int i=0;i<Q;i++)
                nnew[ix][i] = n[ix][(i+1)%2];
    }

}


void LatticeGas::Adveccione(void){

    for(int ix=0;ix<Lx;ix++) //Para cada celda
            for(int i=0;i<Q;i++) //Para cada dirección
                n[(ix+Lx+V[i])%Lx][i] = nnew[ix][i];
}



double LatticeGas::rho(int ix,bool UseNew){
  if(UseNew)
    return nnew[ix][0]+nnew[ix][1];
  else
    return n[ix][0]+n[ix][1];
}


double LatticeGas::GetSigma2(void){
    int ix;

    //Calcular cuántas bolitas hay.
    double N=0;

    for(ix=0;ix<Lx;ix++)
        N+=rho(ix,false);

    double xprom = 0;
    for (ix=0;ix<Lx ;ix++ ) {
        xprom += ix*rho(ix,false);
    }
    xprom/=N;


    //Varianza promedio
    double Sigma2=0;

    for(ix=0;ix<Lx;ix++)
        Sigma2+=pow(ix-xprom,2.0)*rho(ix,false);
    Sigma2/=(N-1);

    return Sigma2;
}



int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);

  int N = (int) Lx*0.4;
  double mu=Lx/2;
  double sigma = Lx/8;

  int t, tmax=400;
  Difusion.Borrese();
  Difusion.Inicie(N,mu,sigma,ran64);

  for(t=0;t<tmax;t++){

      cout<<t<<" "<<Difusion.GetSigma2()<<endl;
      Difusion.Colisione(ran64);
      Difusion.Adveccione();
  }



  return 0;
}
