#include <iostream>
#include <math.h>

// n debe ser par.

//(h/3)* (f(x0) + 4*f(x1) + 2*f(x2) + ... + 4(f(xn-1)) + f(xn))


double f(double x){

  return std::cos(x);
}

double simpson(double a, double b, int n){

  n*=2;
  double h= (b-a)/n;
  double suma = 0.0;
  for (int i=0;i<n+1; i++ ) {
    double x =a+i*h;

    if (i==0 || i==n) {
      suma = suma + 2*f(x);
    }
    else if (i%2==0) {
      suma= suma + 2*f(x);
    }
    else {
      suma = suma + 4*f(x);

    }

    return suma*h/3.0;

  }

}

int main(int argc, char *argv[]) {

  std::cout<<"La integral por Simpson es: "<< simpson(0,M_PI/2,14);


        return 0;
}
