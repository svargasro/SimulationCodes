#include <iostream>
#include <math.h>


double f(double t, double x){
  return x/2.0;}

void unPasoRungeKutta4(double &t, double &x, double dt ){

  double dx1, dx2, dx3, dx4;

  dx1 = dt*f(t,x);
  dx2 = dt*f(t+ (dt/2), (x + dx1/2));
  dx3 = dt*f(t+(dt/2),x+(dx2/2));
  dx4= dt*f(t+dt,x+dx3);
  x+=(dx1+2*(dx2+dx3)+dx4)/6.0;
  t+=dt;
}



int main(int argc, char *argv[]) {


  double t,x;
  double dt=0.01;
  for(t=0,x=1;t<10;){

unPasoRungeKutta4(t,x,dt);

  }

 std::cout<<"x(10) = "<<x<<std::endl;
 std::cout<<std::exp(t/2)<<std::endl;

  return 0;
}
