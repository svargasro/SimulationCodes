#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using std::vector;
using std::pow;
using std::sqrt;
using std::ofstream;
using std::endl;
using std::cout;

double GM = 1.0;
double Deltat = 0.001;
int Nsteps= 200000;


class Planeta{
private:
  vector <double> r;
  vector <double> V;
  vector <double> F;
  double m;

public:

  Planeta(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double m0){
  r = {x0, y0, z0};
  V = {Vx0,Vy0,Vz0};
  m=m0;
  F = {0.0,0.0,0.0};
  }



  void CalculeFuerza(void){
    double norma = sqrt(pow(r[0],2) + pow(r[1],2)+ pow(r[2],2));
    double factor = -GM*m/(pow(norma,3));

    for(int i=0; i<=2 ; i++ ){
      F[i] = factor*r[i];
    }



  }

  void Muevase(double dt){
    for (int i=0;i<=2 ; i++) {
    r[i] = r[i] + dt*V[i];

    V[i] = V[i] + (dt/m)*F[i];

    }

  }

  double getX(void){
    return r[0];
  }

  double getY(void){
    return r[1];
  }



};

int main(int argc, char *argv[]) {
  double r0 = 10.0;
  double omega= sqrt(GM/(pow(r0,3)));
  double V0 = omega*r0;

  vector<double> xdata(Nsteps);
  vector<double> ydata(Nsteps);
//(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double m0)
  Planeta planeta1(r0,0,0,0,0.5*V0,0,0.453);





for (int i=0;i<=Nsteps;i++) {

  cout<<planeta1.getX()<<" "<<planeta1.getY()<<endl;
  planeta1.CalculeFuerza();
  planeta1.Muevase(Deltat);

 }

  return 0;
}
