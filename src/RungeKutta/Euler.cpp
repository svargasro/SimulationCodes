#include <iostream>
#include <math.h>


double f(double t, double x){

  return (x/2.0);
}

void unPasoEuler(double &t, double &x, double dt)
{
  double dx;
  dx = dt*f(t,x);
  x+= dx;
  t+= dt;
}


int main(int argc, char *argv[]) {


  double t,x;
  double dt=0.01;
  for(t=0,x=1;t<10;){

                  unPasoEuler(t,x,dt);

  }

 std::cout<<"x(10) = "<<x<<std::endl;
 std::cout<<std::exp(t/2)<<std::endl;

               return 0;
}
