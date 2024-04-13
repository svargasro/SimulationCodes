#include <iostream>
#include <cmath>
#include "vector.h"
#include <fstream>
using namespace std;

//Constantes globales

const int N=5;
const double g=980;
const double L = 12.0;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
	double Theta, Omega,Tau;  double m,R, I, L, x0;
public:
  void Inicie(double Theta0,double Omega0,double m0, double R0, double L0, double x00);
  void BorreTorque(void){Tau=0;};
  void SumeTorque(double Tau0){Tau+=Tau0;};
  void Mueva_Theta(double dt,double coeficiente);
  void Mueva_Omega(double dt,double coeficiente);
  void Dibujese(void);
	double GetX(void){return x0 +L*sin(Theta);};
	double GetY(void){return -L*cos(Theta);};
	double GetTheta(void){return Theta;};
  double GetTau(void){return Tau;}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double Theta0,double Omega0,double m0, double R0, double L0, double x00){
  Theta = Theta0; Omega = Omega0; m=m0; R= R0; L= L0; I =m*L*L; x0 = x00;
}
void Cuerpo::Mueva_Theta(double dt,double coeficiente){
  Theta += Omega*( coeficiente*dt);
}
void Cuerpo::Mueva_Omega(double dt,double coeficiente){
  Omega += Tau*( coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<GetX()<<"+"<<R<<"*cos(t),"<<GetY()<<"+"<<R<<"*sin(t)";
	 cout<<" , "<<x0<<"+"<<L/7<<"*t*sin("<<Theta<<"),-"<<L/7<<"*t*cos("<<Theta<<")";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Pendulo, double K);
  void CalculeFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2, double K);    
};
void Colisionador::CalculeFuerzas(Cuerpo * Pendulo, double K){
  int i,j;
  //Borrar todas las fuerzas
  for(i=0;i<N;i++)
    Pendulo[i].BorreTorque();
   //Adicionar fuerzas individuales (que no son de interaccion)
   for(i=0;i<N;i++)
    Pendulo[i].SumeTorque(-Pendulo[i].L*Pendulo[i].m*g*sin(Pendulo[i].Theta));
  //Calcular las fuerzas entre todas las parejas de Pendulos
  for(i=0;i<N-1;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Pendulo[i],Pendulo[j], K);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2, double K){
    double S = (Pendulo1.GetX() + Pendulo1.R) -(Pendulo2.GetX() - Pendulo2.R);
	if(S>0){
	   double F = K* pow(S, 1.5);
	   double T2= F*Pendulo2.L;
	   Pendulo1.SumeTorque(-T2); Pendulo2.SumeTorque(T2);
	}
}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<"\n"; 
  cout<<"set output 'CunadeNewton.gif'"<<"\n";
  cout<<"unset key"<<"\n";
  cout<<"set xrange[-14:34]"<<"\n";
  cout<<"set yrange[-18:0]"<<"\n";
  cout<<"set size ratio -1"<<"\n";
  cout<<"set parametric"<<"\n";
  cout<<"set trange [0:7]"<<"\n";
  cout<<"set isosamples 12"<<"\n";  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<"\n";
}


double funcion(double x,double x1,double tau1,double x2, double tau2)
{
    double pendiente=((tau2-tau1)/(x2-x1));
	return pendiente*(x-x1)+tau1;
}

double Cerosrecta(double p1x,double p1y,double p2x,double p2y) 
{   
    double a=p1x, b=p2x;
	double  m = 0, fa, fm;
	const double Epsilon = 1e-12;

	fa = funcion(a,p1x,p1y,p2x,p2y);

	while ((b - a) > Epsilon)
	{
		m = (a + b) / 2.0;
		fm = funcion(m,p1x,p1y,p2x,p2y);
		if (fa * fm > 0)
		{
			a = m;
			fa = funcion(a,p1x,p1y,p2x,p2y);
		}
		else
			b = m;
	}
	return (a + b) / 2;
}

int main(){

    double Kc[7]= {0.1e9, 0.2e9, 0.5e9, 1e9, 2e9, 5e9, 10e9};
    string name[7] = {"01","02","05","1","2","5","10"};
    //sistema de maximos
    ofstream maxs;
    maxs.open("periodos.csv");
    double T1,T2,T3,t1,t2,t3;
    //----------
    for(int ii = 0; ii <7; ii++){
    double K = Kc[ii];
    ofstream data;
    data.open ("C"+ name[ii] +".csv");
		
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=100, R0=2, L0=12;
  double T=2*M_PI*sqrt(L0/g);
  double t,tdibujo,tcuadro=T/50;
  double tmax=0.2,dt=1e-5, tstart;
  int i,plot,max,q;

  
  //---------------(Theta0,Omega0,m0,R0,L0,x00)
  Pendulo[0].Inicie(-M_PI/12,0,m0,R0,L0,R0);
	for(i=1; i <N; i++)
			Pendulo[i].Inicie(0,0,m0,R0,L0,(2*i+1)*R0);
  
  T1=Pendulo[1].GetTau();
  t1=0;

  
  InicieAnimacion();
  
  for(t=0,tdibujo=0,plot=0,max=0,q=0; t<tmax; t+=dt,tdibujo+=dt){
    //Dibujar
  if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) 				Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    if(plot==0){
        if(Pendulo[1].GetTau()>300){ 
            if(t>0.1)  {
              tstart=t;
            plot=1;}};}
    if(plot==1)
        {
        if(signbit(T3)!=signbit(T1))
            {   if(q%2==0){
                q++;
                //maxs<<K<<";"<<t3<<";"<<T3<<";"<<t1<<";"<<T1<<"\n";
                maxs<<K<<";"<<Cerosrecta(t3,T3,t1,T1)-tstart<<"\n";
                }
            }
        
        data<<t<<";"<<Pendulo[1].GetTau() <<"\n";
        }
    //
    T3=T2;                  t3=t2;
    T2=T1;                  t2=t1;
    T1=Pendulo[1].GetTau(); t1=t;
    // Mover por PEFRL
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Zeta);
    Newton.CalculeFuerzas(Pendulo,  K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Chi);
    Newton.CalculeFuerzas(Pendulo,  K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Coeficiente2);
    Newton.CalculeFuerzas(Pendulo, K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Chi);
    Newton.CalculeFuerzas(Pendulo, K);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_Theta(dt,Zeta);   
  }
		data.close();
		}
  
  return 0;
}
