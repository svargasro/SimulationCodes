#include <iostream>
#include <cmath>
#include "vector.h"
using std::cout;
using std::endl;



int main(int argc, char *argv[])  {
vector3D a,b,c,d;
 a.load(1,0,0);
 b.load(0,1,0);
 c.load(0,0,1);
 //c=a+b;

//Hace primero la suma, luego la operación lógica
 d= (a^b) + c;
 d.show();

  return 0;
}
