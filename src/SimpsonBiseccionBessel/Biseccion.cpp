#include <iostream>
#include <math.h>
using std::cout;
using std::endl;


const double ERR = 1e-6;

double f(double x){
  return std::sin(x)/x;
}


double CeroPorBiseccion(float a, float b){

  double m,fa,fm;
  //Calcule el cero por bisección de f(x) entre a y b.

  fa=f(a);

  while((b-a)>ERR)  {
    m=(a+b)/2.0;
    fm = f(m);
    if (fm*fa<0) {
      b=m;
    }
    else {
      a=m;
      fa=fm;
    }

  }



  return (a+b)/2;
}

int main(int argc, char *argv[]) {

cout<<"El cero es x="<<CeroPorBiseccion(2,4)<<endl;





  return 0;
}
