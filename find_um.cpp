#include <arprec/mp_real.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "integrate_backwards.cpp"

mp_real H(mp_real x,mp_real y);

using std::cout;
//using std::cerr;
using std::endl;
using namespace std; 


int main(int argc, char *argv[]){
/*
This program just reads some data, calls the integrate routine and prints the resulting value.
In this case, it calls integrate_backwards in order to approach the unstable manifold of \tL^-.
*/

mp::mp_init(100);

mp_real *z0=new mp_real[5];
mp_real *z1=new mp_real[5];
mp_real delta,tau,yh,step;
mp_real stop_tol;
mp_real *error=new mp_real;

double yhmax;
int numiterates;
ofstream fout;

z0[0]=mp_real(argv[1]);
z0[1]=mp_real(argv[2]);
z0[2]=mp_real(argv[3]);
z0[3]=mp_real(argv[4]);
z0[4]=mp_real(argv[5]);


delta=mp_real(argv[6]);


istringstream ( argv[8]) >> numiterates;//Number of iterates to be calculated

stop_tol=mp_real("4");


integrate(z0,z1,delta,numiterates,error);


/*
Compte aquí, que encara que la precissió de sortida sigui tant gran com la de càlcul, hi ha un error
de l'ordre de 1e-20 en imprimir els números.
*/
fout.open(argv[7],ios::out);
//I added this to avoid escaping by tangencies with x=0.
//This tolerance could be up to the order of the error in the Newton method.
if (abs(z1[0])<stop_tol && abs(z1[1]) > mp_real("1e-20")){ 
//if (abs(z1[0])<stop_tol){
	//fout <<0<<" "<<z0[1]<<" "<<z1[0]<<" "<<z1[1]<<" "<<z1[2]<<" "<<z1[3]<<" "<<z1[4]<<" "<<endl;
	fout <<setprecision(30)<<0<<" "<<dble(z0[1])<<" "<<dble(z1[0])<<" "<<dble(z1[1])<<" "<<dble(z1[2])<<" "<<dble(z1[3])<<" "<<dble(z1[4])<<" "<<dble(*error)<<endl;
}
else {
	//fout <<0<<" "<<z0[1]<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<endl;
	//fout <<setprecision(30)<<0<<" "<<dble(z0[1])<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<endl;
	fout <<setprecision(30)<<0<<" "<<dble(z0[1])<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<"Nan"<<" "<<dble(*error)<<endl;
}


fout.close();
delete [] z0;
delete [] z1;
delete error;

return(0);

}

mp_real H(mp_real x,mp_real y){
mp_real h;

h=pow(y,2)/mp_real("2.0")-pow(x,2)/mp_real("2.0");

if (x>0){
	h=h+x;
}
else{
	h=h-x;
}
return h;
}
