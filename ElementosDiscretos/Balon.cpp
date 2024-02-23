#include <iostream>
#include <cmath>
using std::cout;
using std::endl;

const double g=9.8;
class Objeto;


//Implementación clase
class Objeto{
  // private: Es equivalente a tenerlo privado. Es por defecto.
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
  //Declaración funciones.
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);

  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){
    return x;
  }; //PONER ; al final.
  double Gety(void){
    return y;
  }; //Inline corresponde a una macro (Mucho más rápido, porque lo pasa a código.). Por fuera corresponde a una subrutina.

};



//Implementación funciones
void Objeto::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0;
  y=y0;
  Vx=Vx0;
  m=m0;
  Vy=Vy0;
  R=R0;
} //Subrutina (Compilado más corto pero demora un poco más en ejecutar).

void Objeto::CalculeFuerza(void){
  Fx=0;
  Fy=-m*g;
}

void Objeto::Muevase(double dt){
  x+=Vx*dt;
  y+=Vy*dt;
  Vx+= Fx*dt/m;
  Vy+= Fy*dt/m;

}

void Objeto::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
  //cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<theta<<"),-"<<l/7<<"*t*cos("<<theta<<")";
}

//Funciones globales

//---FuncionesDeAnimacion


void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1:41]"<<endl;
  cout<<"set yrange[-10:8]"<<endl;
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





int main(int argc, char *argv[]){

  double t,dt=0.001;
  Objeto Balon;
  InicieAnimacion();
  Balon.Inicie(0,0,16,9,0.453,2);

  for (t=0;t<2.5;t+=dt) {
    InicieCuadro();
    Balon.Dibujese();
    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
                                           Balon.CalculeFuerza();
                                           Balon.Muevase(dt);
                                           TermineCuadro();
  }
                        return 0;
}
