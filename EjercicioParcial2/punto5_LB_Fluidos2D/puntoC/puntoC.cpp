#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

const int Lx=512;
const int Ly=64;

const int Q=9;

const double tau=1.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  double w[Q];      //Weights 
  int Vx[Q],Vy[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Collision(void);
  void ImposeFields(double Ufan, int ixc, int iyc,int R);
  void Advection(void);
  vector<double> Derivatives(int ix, int iy, double dt); //Calcula dUx/dx, dUx/dy, dUy/dx, dUy/dy dada una celda. 

  //Ya con el valor de las derivadas, se calculan las componentes del tensor de esfuerzos.
  double sigmaxx(double rho0, double eta, double dx_dx); 
  double sigmaxy(double rho0, double eta, double dx_dy, double dy_dx);
  double sigmayy(double rho0, double eta, double dy_dy);

  //Función de interpolación bilineal
  vector<double> interpolationSigma(double x, double y, double nu, double dt);

  //Función que calcula el diferencial de fuerza
  vector<double> dF(double x, double y, double dx, double dy, double nu, double dt);

  
  void Start(double rho0,double Ux0,double Uy0);
  void Print(const char * NameFile,double Ufan);
};  
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;

              Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
              Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Start(double rho0,double Ux0,double Uy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ix,iy,i);
	f[n0]=feq(rho0,Ux0,Uy0,i);
      }
}  
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){ //for each velocity vector
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }  
}

void LatticeBoltzmann::ImposeFields(double Ufan, int ixc, int iyc,int R){
  int i,ix,iy,n0; double rho0; double R2=R*R;
  //go through all cells, looking if they are fan or obstacle
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //fan
      if(ix==0)
	for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ufan,0,i);}
      //obstacle
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2) 
	for(i=0;i<Q;i++) {n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}
      //An extra point at one side to break the isotropy
      // else if(ix==ixc && iy==iyc+R+1)
      //for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}	
    }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
vector<double> LatticeBoltzmann::Derivatives(int ix, int iy,double dt){
   int i,ixnext,iynext;
   double rhoNext,UxNext,UyNext;
   double sumxx=0, sumxy=0, sumyx=0, sumyy=0;

  for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;

  //Calculo las velocidades macroscópicas para las celdas aledañas a la celda en la que estoy parada

	rhoNext=rho(ixnext,iynext,false); UxNext=Jx(ixnext,iynext,false)/rhoNext; UyNext=Jy(ixnext,iynext,false)/rhoNext; //preguntarse si usar false o true: false si se usa luego de advección.

  //Posibles combinaciones de la sumatoria para calcular más adelante las derivadas (se muestra en el taller)

	sumxx += w[i]*Vx[i]*UxNext;
	sumxy += w[i]*Vy[i]*UxNext; 
	sumyx += w[i]*Vx[i]*UyNext;
	sumyy += w[i]*Vy[i]*UyNext;
	
  }

  //Ya las posibles derivadas calculadas

	double dx_dx = (3.0/dt)*sumxx;
  double dx_dy = (3.0/dt)*sumxy; //dx_dy es distinto a dy_dx 
	double dy_dx = (3.0/dt)*sumyx;
	double dy_dy = (3.0/dt)*sumyy;

  //Vector que guarda las derivadas

  vector<double> derivs = {dx_dx, dx_dy, dy_dx, dy_dy};
  return derivs;
}

//Calcula la componente xx del tensor de esfuerzos
double LatticeBoltzmann::sigmaxx(double rho0, double eta, double dx_dx){
  double p = rho0/3.0;
  double sxx = -p+eta*(2*dx_dx);
  return sxx;
}

//Calcula la componente xy e yx (son iguales) del tensor de esfuerzos
double LatticeBoltzmann::sigmaxy(double rho0, double eta, double dx_dy, double dy_dx){
  double sxy = eta*(dx_dy + dy_dx);
  return sxy;
}
//Calcula la componente yy del tensor de esfuerzos
double LatticeBoltzmann::sigmayy(double rho0, double eta, double dy_dy){
  double p = rho0/3.0;
  double syy = -p+eta*(2*dy_dy);
  return syy;
}

vector<double> LatticeBoltzmann::interpolationSigma(double x, double y, double nu, double dt){
/*
M*ancho = x + 0  ;   x: punto cualquiera que mido en el plano 
M = x/ancho

Lx=5
0<=ix<=4

|__0__|__1__|__2__|__3__|__4__|

x=2.5, ¿A qué ix corresponde?
M=ix= int(2.5/ancho) = 2
*/

int ix = static_cast<int>(x); 
int iy = static_cast<int>(y); 
double u = x - ix;
double v = y - iy;

//Cálculo para celda (ix,iy)
double rho_ix_iy = rho(ix,iy,false);
double eta_ix_iy = nu*rho_ix_iy;

vector<double> derivatesIxIy = Derivatives(ix, iy, dt); //{dx_dx, dx_dy, dy_dx, dy_dy}
double sigmaxx_ix_iy = sigmaxx(rho_ix_iy, eta_ix_iy, derivatesIxIy[0]);
double sigmaxy_ix_iy = sigmaxy(rho_ix_iy, eta_ix_iy, derivatesIxIy[1],derivatesIxIy[2]);
double sigmayy_ix_iy = sigmayy(rho_ix_iy, eta_ix_iy, derivatesIxIy[3]);

//Cálculo para (ix+1,iy)
double rho_ixP1_iy = rho(ix+1,iy,false);
double eta_ixP1_iy = nu*rho_ixP1_iy;

vector<double> derivatesIxP1Iy = Derivatives(ix+1, iy, dt); //{dx_dx, dx_dy, dy_dx, dy_dy}
double sigmaxx_ixP1_iy = sigmaxx(rho_ixP1_iy, eta_ixP1_iy, derivatesIxP1Iy[0]);
double sigmaxy_ixP1_iy = sigmaxy(rho_ixP1_iy, eta_ixP1_iy, derivatesIxP1Iy[1],derivatesIxP1Iy[2]);
double sigmayy_ixP1_iy = sigmayy(rho_ixP1_iy, eta_ixP1_iy, derivatesIxP1Iy[3]);

//Cálculo para (ix,iy+1) [Problema si iy es el límite.]
double rho_ix_iyP1 = rho(ix,iy+1,false);
double eta_ix_iyP1 = nu*rho_ix_iyP1;

vector<double> derivatesIxIyP1 = Derivatives(ix, iy+1, dt); //{dx_dx, dx_dy, dy_dx, dy_dy}
double sigmaxx_ix_iyP1 = sigmaxx(rho_ix_iyP1, eta_ix_iyP1, derivatesIxIyP1[0]);
double sigmaxy_ix_iyP1 = sigmaxy(rho_ix_iyP1, eta_ix_iyP1, derivatesIxIyP1[1],derivatesIxIyP1[2]);
double sigmayy_ix_iyP1 = sigmayy(rho_ix_iyP1, eta_ix_iyP1, derivatesIxIyP1[3]);


//Cálculo para (ix+1,iy+1) [Problema si iy es el límite.]
double rho_ixP1_iyP1 = rho(ix+1,iy+1,false);
double eta_ixP1_iyP1 = nu*rho_ixP1_iyP1;

vector<double> derivatesIxP1IyP1 = Derivatives(ix+1, iy+1, dt); //{dx_dx, dx_dy, dy_dx, dy_dy}
double sigmaxx_ixP1_iyP1 = sigmaxx(rho_ixP1_iyP1, eta_ixP1_iyP1, derivatesIxP1IyP1[0]);
double sigmaxy_ixP1_iyP1 = sigmaxy(rho_ixP1_iyP1, eta_ixP1_iyP1, derivatesIxP1IyP1[1],derivatesIxP1IyP1[2]);
double sigmayy_ixP1_iyP1 = sigmayy(rho_ixP1_iyP1, eta_ixP1_iyP1, derivatesIxP1IyP1[3]);


double interpolatedSigmaXX = sigmaxx_ix_iy*(1-u)*(1-v) + sigmaxx_ixP1_iy*u*(1-v) + sigmaxx_ix_iyP1*(1-u)*v + sigmaxx_ixP1_iyP1*u*v;
double interpolatedSigmaXY = sigmaxy_ix_iy*(1-u)*(1-v) + sigmaxy_ixP1_iy*u*(1-v) + sigmaxy_ix_iyP1*(1-u)*v + sigmaxy_ixP1_iyP1*u*v;
double interpolatedSigmaYY = sigmayy_ix_iy*(1-u)*(1-v) + sigmayy_ixP1_iy*u*(1-v) + sigmayy_ix_iyP1*(1-u)*v + sigmayy_ixP1_iyP1*u*v;

vector <double> interpolatedSigma = {interpolatedSigmaXX, interpolatedSigmaXY, interpolatedSigmaYY};
return interpolatedSigma;

 
}
vector<double> LatticeBoltzmann::dF(double x,double y, double dx, double dy, double nu, double dt){


vector<double> sigmaInterpolated = interpolationSigma(x,y,nu,dt);

vector<double> fuerza = {sigmaInterpolated[0]*dx+sigmaInterpolated[1]*dy,sigmaInterpolated[1]*dx+sigmaInterpolated[2]*dy};

return fuerza; 

/*
dFi = sigma(ij) dAj
-dFx = sigma(xx)dAx + sigma(xy)dAy
-dFy = sigma(yx)dAx + sigma(yy)dAy

dA = (dx,dy)
*/

}




void LatticeBoltzmann::Print(const char * NameFile,double Ufan){
  ofstream MyFile(NameFile); int ix,iy; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix+=4){
    for(iy=0;iy<Ly;iy+=4){
      rho0 = rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
      MyFile<<ix<<" "<<iy<<" "<<Ux0/Ufan*4<<" "<<Uy0/Ufan*4<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//--------------- Global Functions ------------

int main(void){
  LatticeBoltzmann Air;
  int t,tmax=5000;
  double rho0=1.0,Ufan0=0.1;
  double dt = 0.1;
  double nu = dt*(1/3.0)*(tau- 1.0/2);
  int ixc=128, iyc=32, R=8;

  //Start
  Air.Start(rho0,Ufan0,0);
  //Run
  for(t=0;t<tmax;t++){
    Air.Collision();
    Air.ImposeFields(Ufan0,ixc,iyc,R);
    Air.Advection();
  }
  //Show
  Air.Print("wind.dat",Ufan0);
 
  return 0;
}  
