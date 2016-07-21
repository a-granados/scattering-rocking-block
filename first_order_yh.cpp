#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
/*
We calculate here the intersection of the heteroclinic connection with x=0 up to third order in epsilon:

y^h_eps=y^h_0+eps*\int_{-\infty}^0{X^-,h}(unperturbed_heteroclinic)dt+O(eps^2)
This value depends on eps, tau,v and s.
The obtained value could be use in find_het.sh to accelerate to process, or to check that everything worked fine.
*/

using std::cout;
//using std::cerr;
using std::endl;
using namespace std;

double M(double gamm,double tau0,double v0,double s0,double params[]);
double PoisBra(double t, double gamm,double tau0,double v0, double s0,double params[]);
double rmod(double x,double y);

int main (int argc, char *argv[]){
cout.precision(30);

double tau0,s0,v0,h,Mval,ystar,y0h;
double delta;
double params[3];
ifstream fin;

fin.open("system_params.dat",ios::in);
fin>> params[0];//omega
fin>> params[1];//te
fin>> params[2];//tk
fin.close();


tau0=atof(argv[1]);
v0=atof(argv[2]);
s0=atof(argv[3]);
delta=atof(argv[4]);
y0h=1;

Mval=M(0.0,tau0,v0,s0,params);
ystar=y0h+delta*params[1]/y0h*Mval;

cout<<ystar<<endl;

}

double M(double gamm,double tau0,double v0,double s0,double params[]){
double sum,h;
int i,imax;

h=1e-4;
imax=100000;


sum=0;

for (i=1;i<imax;i++){
	//Simpson integration
//Unlike in the Melnikov function, we integrate here only backwards.
	sum=sum+PoisBra(i*h,gamm,tau0,v0,s0,params)*h;
//	sum=sum+PoisBra(-i*h,gamm,tau0,v0,s0,params)*h;

}

return sum;
}

double PoisBra(double t, double gamm,double tau0,double v0, double s0,double params[]){

double alpha,alphap,ft;
double x,y,u,v,s,A;
double omega,te,tk;

omega=params[0];
te=params[1];
tk=params[2];

alphap=log((1+v0)/(1-v0));
alpha=2*alphap;

A=t+tau0+gamm;

A=rmod(A,alpha);
if (A<alphap){
	u=(v0-1)/2*exp(A)-(v0+1)/2*exp(-A)+1;
	v=(v0-1)/2*exp(A)+(v0+1)/2*exp(-A);
}
else{
	u=(-v0+1)/2*exp(-alphap+A)+(v0+1)/2*exp(alphap-A)-1;
	v=(-v0+1)/2*exp(-alphap+A)-(v0+1)/2*exp(alphap-A);
}

if (t>=0){
	x=1-exp(-t);
	y=exp(-t);
}
else{
	x=exp(t)-1;
	y=exp(t);
}

s=s0+t+gamm;
//ft=sin(omega*s)+pow(sin(omega*s),2);
ft=cos(omega*s);


//cout <<omega<<endl;
//cout <<x<<" "<<y<<" "<<u<<" "<<v<<" "<<s<<endl;

//return -y*(te*ft+tk*(x-u));
 //Perturbation which vanishes at Lambda 
return -4*y*x*(x-1)*(x+1)*u*ft;

}

double rmod(double x,double y){
//Return x mod y, with x,y>0 real numbers.

int i;

if (x>=0){
i=1;
while (i*y<x){
	i++;
}
return x-(i-1)*y;
}
else {
i=-1;
while (x<i*y){
	i=i-1;
}
//return (i+1)*y-x;
return x-i*y;
}
}


