#include <iostream>
#include <math.h>


const double T =M_PI;
const double omega = 2*M_PI/T;
const double omega2=omega*omega;

double f1(double t, double x1, double x2){
  return -omega2*x2;}
double f2(double t, double x1, double x2){
return x1;
}

void unPasoRungeKutta4(double &t, double &x1,double &x2, double dt ){

  double dx11, dx21, dx31, dx41;
  double dx12, dx22, dx32, dx42;

  dx11 = dt*f1(t,x1,x2);
  dx12 = dt*f2(t,x1,x2);

  dx21 = dt*f1(t+ (dt/2), (x1 + dx11/2), x2+dx12/2);
  dx22 = dt*f2(t+ (dt/2), (x1 + dx11/2), x2+dx12/2);


  dx31 = dt*f1(t+(dt/2),x1+(dx21/2), x2+ (dx22/2));
  dx32 = dt*f2(t+(dt/2),x1+(dx21/2), x2+ (dx22/2));

  dx41 = dt*f1(t+dt,x1+dx31, x2+dx32);
  dx42 = dt*f2(t+dt,x1+dx31, x2+dx32);


  x1+=(dx11+2*(dx21+dx31)+dx41)/6.0;
  x2+=(dx12+2*(dx22+dx32)+dx42)/6.0;
  t+=dt;
}



int main(int argc, char *argv[]) {


  double t,x1,x2;
  double dt=0.01;
  for(t=0,x1=1, x2=0;t<10;){

    unPasoRungeKutta4(t,x1,x2,dt);
    std::cout<<t<<" "<<x2<<std::endl;

  }




  return 0;
}
