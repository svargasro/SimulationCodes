#include <iostream>
#include <math.h>

double f(double t, double x, int n){

  return std::cos(n*t - x*std::sin(t));
}

double IntegralPorSimpson(double a, double b, int nsteps, double x, int n){

  nsteps*=2;
  double h=(b-a)/nsteps;
  double suma=0;

  for (int i=0 ; i<(nsteps+1) ; i++ ) {
    double t=a+i*h;

    if (i==0 || i==nsteps) {
      suma+= f(t,x,n);
    }
    else if (i%2 == 0) {
      suma+=2*f(t,x,n);
    }
    else{
      suma+=4*f(t,x,n);
    }

  }

  return suma*h/3;

}

double Bessel(int n, double x){
  return (1/M_PI)*IntegralPorSimpson(0,M_PI,100,x,n);
}





int main(int argc, char *argv[]) {
  int n =3;

  //int n= std::stoi(argv[1]);

  //DatosParaGraficar
  double h = 100;
  double a = 0;
  double b = 10;
  double dx = (b-a)/h;

  for (double x=a;x<b ; x=x+dx) {

    std::cout<<x<<" "<<Bessel(n,x)<<std::endl;
  }



             return 0;
}
