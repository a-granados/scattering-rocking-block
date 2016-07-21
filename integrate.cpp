#include <arprec/mp_real.h>
#include <iostream>
#include <fstream>
#include <sstream>

void integrate(mp_real *xin,mp_real *xout,mp_real delta,int poincmapiters,mp_real *error);
void f(mp_real *x, mp_real *dx, mp_real params[], mp_real delta, int lblock, int rblock);
void runge4(mp_real *xprev, mp_real *xnext, mp_real params[],mp_real delta, int lblock, int rblock, mp_real step);
void rk78(mp_real *x,mp_real *xnext, mp_real params[],mp_real delta, int lblock, int rblock,mp_real *ah, mp_real tol,mp_real *totalerror);
int get_next_impact(mp_real *pimpact,mp_real *nimpact,mp_real *accumTime,mp_real params[],mp_real delta,int *lblock, int *rblock,mp_real *error);
mp_real SGN(mp_real x);
mp_real MAX(mp_real a,mp_real b);
void obtain_constants(mp_real *alfa,mp_real *beta,mp_real *c7,mp_real *c8);

/*
This is a modification of the program "integrate" in "fixed_step_size_RK4" folder.
*/
//#include <string>

using std::cout;
//using std::cerr;
using std::endl;
using namespace std; 

#define N 5
#define plotorbit true


void integrate(mp_real *xin,mp_real *xout,mp_real delta,int poincmapiters, mp_real *error){
/*
Call:
integrate(zin,zout,delta,num_stroboscopic_iterates_to_simulate)
The value after poincmapiters iterations of the poincare map is stored in x.
*/

mp::mp_init(100);
cout.precision(30);

mp_real *x=new mp_real[5];
mp_real *y=new mp_real[5];
mp_real params[4],*accumTime;
int i,*lblock,*rblock,flag;
int lblockini;
double tf;
ifstream fin;
ofstream fout;
//const char* tmp;
string tmp;

x[0]=xin[0];
x[1]=xin[1];
x[2]=xin[2];
x[3]=xin[3];
x[4]=xin[4];

lblock=new int;
rblock= new int;
accumTime=new mp_real;
mp_real pi = mp_real::_pi;


fin.open("system_params.dat",ios::in);
fin>> tmp;
params[0]=mp_real(tmp);//omega
fin>> tmp;
params[1]=mp_real(tmp);//te
fin>> tmp;
params[2]=mp_real(tmp);//tk
fin.close();

if (x[0]>0){
	*lblock=1;
}
else if (x[0]<0){
	*lblock=-1;
}
else{
	if (x[1]>0){
	*lblock=1;
	}
	else{
	*lblock=-1;
	}
}

lblockini=*lblock; //This is to check whether we cross x=0 or not

if (x[2]>0){
	*rblock=1;
}
else if (x[2]<0){
	*rblock=-1;
}
else {
	if (x[3]>0){
	*rblock=1;
	}
	else{
	*rblock=-1;
	}
}
/////////////////////
/////////////////////

/*
mp_real *ah;
ah=new mp_real;
*ah=mp_real("1e-3");
cout<<"hello"<<endl;
rk78(x,y,params,mp_real("0.05"),1,1,ah,mp_real("1e-5"));
delete ah;
cout<<"hello"<<endl;
exit(0);
*/
*accumTime=mp_real("0");
//fout.open(argv[7],ios::out);

//
///And we start the iterative process:
//

i=0;
*error=0;
//fout<<dble(x[0])<<" "<<dble(x[1])<<" "<<dble(x[2])<<" "<<dble(x[3])<<" "<<dble(x[4])<<" "<<i<<endl;///We save the initial condition


if (plotorbit) {
	cout<<dble(x[0])<<" "<<dble(x[1])<<" "<<dble(x[2])<<" "<<dble(x[3])<<" "<<dble(x[4])<<" "<<i<<endl;
}

flag=get_next_impact(x,y,accumTime,params,delta,lblock,rblock,error);
//We check whether an image of the stroboscopic map has occured.
if (flag==1){
	*accumTime=0;
	i=1;
//	fout<<dble(y[0])<<" "<<dble(y[1])<<" "<<dble(y[2])<<" "<<dble(y[3])<<" "<<dble(y[4])<<" "<<i<<endl;
//	cout <<"Processing the "<<i+1<<"th of "<<poincmapiters<<endl;
}

if (plotorbit){
	if (flag==0){
		cout<<dble(y[0])<<" "<<dble(y[1])<<" "<<dble(y[2])<<" "<<dble(y[3])<<" "<<dble(y[4])<<endl;
	}
}


//To find the heteroclinic connection we also stop when the switching manifold x=0 is crossed.
//while (i<poincmapiters && abs(y[0])<mp_real("2") && abs(y[1])<mp_real("2") && *lblock*lblockini>0){
//while (i<poincmapiters && flag!=2){
while (i<poincmapiters && flag!=2 && *lblock*lblockini>0){
//while (i<poincmapiters && abs(y[1])<mp_real("10")){
	x[0]=y[0];
	x[1]=y[1];
	x[2]=y[2];
	x[3]=y[3];
	x[4]=y[4];
	flag=get_next_impact(x,y,accumTime,params,delta,lblock,rblock,error);
	if (flag==1){
		i=i+1;
//		fout<<dble(y[0])<<" "<<dble(y[1])<<" "<<dble(y[2])<<" "<<dble(y[3])<<" "<<dble(y[4])<<" "<<i<<endl;
		*accumTime=0;
//		cout <<"Processing the "<<i+1<<"th of "<<poincmapiters<<endl;
	}

	if (plotorbit){
		if (flag==0){
			cout<<dble(y[0])<<" "<<dble(y[1])<<" "<<dble(y[2])<<" "<<dble(y[3])<<" "<<dble(y[4])<<endl;
		}
	}

//cout << "Accumulated error: "<<dble(*error)<<endl;
}
	xout[0]=y[0];
	xout[1]=y[1];
	xout[2]=y[2];
	xout[3]=y[3];
	xout[4]=y[4];
//////////////////
//////////////////

//fout.close();
delete [] y;
delete [] x;
delete lblock;
delete rblock;
delete accumTime;

mp::mp_finalize();

}

int get_next_impact(mp_real *pimpact,mp_real *nimpact,mp_real *accumTime,mp_real params[],mp_real delta,int *lblock, int *rblock,mp_real *error){
/*
Given a point pimpact (not necessary at any boundary!) this function returns at nimpact the next impact with the first boundary or the point of the stroboscopic Poincaré map if the crossing is given in the section in time t=0 mod 2*pi/omega.
If the impact with the switching manifold is returned, the value 0 is returned, and 1 otherwise. It returns 2 if the trajectory scapes from the nhim manifold.
This is done by integrating the system, using rung4, until something occurs: one block impacts with the boundary or a zero-crossin in time occurs. In the first case we use a Newton method using the fact that runge4 returns in z[3] the derivative of z[2] at t=z[4]+step (similarly for z[0] and z[1]).
In the variable accumTime we have the accumulated time since the last iterate of the stroboscopic Poincaré map.
*/

mp_real *zaux,step,plotstep;
double tol;
int i, maxiter,block;
mp_real pi = mp_real::_pi;
mp_real stop_tol;
mp_real *ah,rktol;

zaux= new mp_real[N];
ah=new mp_real;
maxiter=2000;
tol=1e-35;//For the zero crossing
step=mp_real("1e-4");
stop_tol=mp_real("4"); //To decide divergence
*ah=mp_real("1e-4");//Initial variable step for the rk78 (it might be varied)
rktol=mp_real("1e-30"); //Maximum absolute error for the rk78

i=1;
//runge4(pimpact,nimpact,params,delta,*lblock,*rblock,step);
rk78(pimpact,nimpact,params,delta,*lblock,*rblock,ah,rktol,error);

/*
if (*accumTime+step>=2*pi/params[0]){//Poincaré stroboscopic map.
//	*accumTime=*accumTime-step;
	step=2*pi/params[0]-(*accumTime);	
	runge4(pimpact,nimpact,params,delta,*lblock,*rblock,step);
	*accumTime=*accumTime+step;
	return 1;
}
else{
*accumTime=*accumTime+step;
}
*/
*accumTime=*accumTime+*ah;


zaux[0]=nimpact[0];
zaux[1]=nimpact[1];
zaux[2]=nimpact[2];
zaux[3]=nimpact[3];
zaux[4]=nimpact[4];

plotstep=step;//This is used to obtain points to plot the trajectories. It only has effect if plotorbit is set to true.
//while (*lblock*zaux[0]>=0 && *rblock*zaux[2]>=0){//We iterate until one rock hits the switching manifold.
//while (*lblock*zaux[0]>=0 && *rblock*zaux[2]>=0 && abs(nimpact[0])<mp_real("2.0")){//We iterate until one rock hits the switching manifold.
//We also stop the integration if we go too away from the hyperbolic point.
//At the moment, the feature the capture the values of the stroboscopic map is
//disabled.
while (*lblock*zaux[0]>=0 && *rblock*zaux[2]>=0 && abs(nimpact[0])<stop_tol){
//while (*rblock*zaux[2]>=0 && abs(nimpact[0])<stop_tol){
	zaux[0]=nimpact[0];
	zaux[1]=nimpact[1];
	zaux[2]=nimpact[2];
	zaux[3]=nimpact[3];
	zaux[4]=nimpact[4];
	//runge4(zaux,nimpact,params,delta,*lblock,*rblock,step);
	rk78(zaux,nimpact,params,delta,*lblock,*rblock,ah,rktol,error);
	if (plotorbit){
		plotstep=plotstep+step;
		if (plotstep>=1e-2){
	cout<< dble(nimpact[0])<<" "<<dble(nimpact[1])<<" "<<dble(nimpact[2])<<" "<<dble(nimpact[3])<<" "<<dble(nimpact[4])<<endl;
		plotstep=0;
		}
	}

	/*
	if (*accumTime+step>=2*pi/params[0]){//Poincaré stroboscopic map.
		step=2*pi/params[0]-(*accumTime);
		runge4(zaux,nimpact,params,delta,*lblock,*rblock,step);
		*accumTime=*accumTime+step;
		return 1;
	}
	else{
	*accumTime=*accumTime+step;
	}
	*/
	*accumTime=*accumTime+*ah;
	i=i+1;
}
if (abs(nimpact[0])>=stop_tol){
	return 2;

}
if (*lblock*zaux[0]<0){//The left block hits first
	block=0;
}
else if (*rblock*zaux[2]<0){//The right block hits first
	block=2;
}
else if (*rblock*zaux[2]<0 && *lblock*zaux[0]<0){
	cout <<"Sa matao Paco"<<endl;
	exit (0);
}

//cout<< "Number of steps to reach the boundary:"<< " " <<i<<endl;

///////////////
///We start now the Newton method
///////////////
i=1;
while (abs(nimpact[block])>tol && i<maxiter){
	step=-nimpact[block]/nimpact[block+1];
	zaux[0]=nimpact[0];
	zaux[1]=nimpact[1];
	zaux[2]=nimpact[2];
	zaux[3]=nimpact[3];
	zaux[4]=nimpact[4];
	runge4(zaux,nimpact,params,delta,*lblock,*rblock,step);
	i=i+1;
/*
	if (*accumTime+step>=2*pi/params[0]){//Poincaré stroboscopic map.
		cout <<"Oju!"<<endl;
		step=2*pi/params[0]-(*accumTime);
		runge4(zaux,nimpact,params,delta,*lblock,*rblock,step);
		*accumTime=*accumTime+step;
		return 1;
	}
	else{
	*accumTime=*accumTime+step;
	}
*/
	*accumTime=*accumTime+step;

}
if (i>=maxiter){
	cout<< "Careful, Newton method in get_next_impact diverges"<<endl;
}
//cout<< dble(nimpact[0])<<" "<<dble(nimpact[1])<<" "<<dble(nimpact[2])<<" "<<dble(nimpact[3])<<" "<<dble(nimpact[4])<<endl;
//cout <<"Number of iterations:"<<" "<< i<<endl;
if (block==0){
	*lblock=-*lblock;
}
else{
	*rblock=-*rblock;
}

delete [] zaux;
delete  ah;
return 0;
}



/*
void get_next_impact(mp_real *pimpact,mp_real *nimpact,mp_real params[],mp_real delta,int lblock, int rblock){
///We try to use here a Bolzano's method
////This is still uncompleted....

mp_real *z1,*z2,*z3,pstep,nstep;
double tol;
int i, maxiter;

z1= new mp_real[N];
z2 = new mp_real[N];
z3 = new mp_real[N];

maxiter=2000;
tol=1e-20;
pstep=mp_real("1e-5");

runge4(pimpact,z2,params,delta,lblock,rblock,pstep);
while (rblock*z2[2]>0){
	z1=z2;;
//	cout << dble(y[0]) << " " << dble(y[1])<<" " << dble(y[2])<<" " << dble(y[3])<<" " << dble(y[4]) << endl;	
	runge4(z1,z2,params,delta,lblock,rblock,pstep);

}

///We start now the Bolzanos method
i=1;

nstep=-pstep/mp_real("2.0");///We know that in the first iteration we have to go backwards.
runge4(z2,z3,params,delta,lblock,rblock,nstep);

while (abs(z3[2])>tol && i<maxiter){
	cout << dble(z3[2])<<endl;
	if (z2[2]*z3[2]>0){
		z1=z2;
		z2=z3;
		pstep=nstep;
		nstep=-pstep/mp_real("2.0");
		runge4(z2,z3,params,delta,lblock,rblock,nstep);
	}
	else {
		z2=z3;
		pstep=nstep;
		nstep=pstep/mp_real("2.0");
		runge4(z2,z3,params,delta,lblock,rblock,nstep);
	}
	i=i+1;
}

delete [] z1;
delete [] z2;
delete [] z3;

}
*/


void f(mp_real *x, mp_real *dx, mp_real params[], mp_real delta, int lblock, int rblock)
{
mp_real omega,te,tk;

omega=params[0];
te=params[1];
tk=params[2];

dx[0]=x[1];
dx[2]=x[3];
dx[4]=1;

if (lblock<0){
//	dx[1]=x[0]+1+delta*(tk*(x[2]-x[0])-te*sin(omega*x[4]));
	dx[1]=x[0]+1+delta*(tk*(x[2]-x[0])-te*cos(omega*x[4]));
//	dx[1]=x[0]+1-delta*4*(x[0]*x[0]-1)*x[0]*x[2]*cos(omega*x[4]);
}
else{
//	dx[1]=x[0]-1+delta*(tk*(x[2]-x[0])-te*sin(omega*x[4]));
	dx[1]=x[0]-1+delta*(tk*(x[2]-x[0])-te*cos(omega*x[4]));
//	dx[1]=x[0]-1-delta*4*(x[0]*x[0]-1)*x[0]*x[2]*cos(omega*x[4]);
}
if (rblock<0){
//	dx[3]=x[2]+1+delta*tk*(x[0]-x[2]);
	dx[3]=x[2]+1+delta*(tk*(x[0]-x[2])-te*cos(omega*x[4]));
//	dx[3]=x[2]+1-delta*(x[0]-1)*(x[0]-1)*(x[0]+1)*(x[0]+1)*cos(omega*x[4]);
}
else{
//	dx[3]=x[2]-1+delta*tk*(x[0]-x[2]);
	dx[3]=x[2]-1+delta*(tk*(x[0]-x[2])-te*cos(omega*x[4]));
//	dx[3]=x[2]-1-delta*(x[0]-1)*(x[0]-1)*(x[0]+1)*(x[0]+1)*cos(omega*x[4]);
}
}



void runge4(mp_real *xprev, mp_real *xnext, mp_real params[],mp_real delta, int lblock, int rblock, mp_real step)
{
//I am assuming that the system is autonomous, and hence t has been included as state variable.


int i;
mp_real *dx1,*dx2,*dx3,*dxprev;
mp_real *x1,*x2,*x3;
mp_real *k1,*k2,*k3,*k4;

dxprev=new mp_real[N];
dx1=new mp_real[N];
dx2=new mp_real[N];
dx3=new mp_real[N];
x1=new mp_real[N];
x2=new mp_real[N];
x3=new mp_real[N];
k1=new mp_real[N];
k2=new mp_real[N];
k3=new mp_real[N];
k4=new mp_real[N];

f(xprev,dxprev,params,delta,lblock,rblock);

for (i=0;i<N;i++)
{
	k1[i]=step*dxprev[i];
	x1[i]=xprev[i]+mp_real("0.5")*k1[i];
}
f(x1,dx1,params,delta,lblock,rblock);
for (i=0;i<N;i++){
	k2[i]=step*dx1[i];
	x2[i]=xprev[i]+mp_real("0.5")*k2[i];
}
f(x2,dx2,params,delta,lblock,rblock);
for (i=0;i<N;i++){
	k3[i]=step*dx2[i];
	x3[i]=xprev[i]+ k3[i];
}
f(x3,dx3,params,delta,lblock,rblock);
for (i=0;i<N;i++){
	k4[i]=step*dx3[i];
}
for (i=0;i<N;i++){
	xnext[i]=xprev[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/mp_real("6.0");
}

delete [] dxprev;
delete [] dx1;
delete [] dx2;
delete [] dx3;
delete [] x1;
delete [] x2;
delete [] x3;
delete [] k1;
delete [] k2;
delete [] k3;
delete [] k4;

}


/*
mp_real rk78(mp_real *at, mp_real * x, mp_real *ah, mp_real tol,
            mp_real hmin, mp_real hmax, int n,
            void (*deriv)(mp_real, mp_real *, int, mp_real *))
*/
void rk78(mp_real *x,mp_real *xnext, mp_real params[],mp_real delta, int lblock, int rblock,mp_real *ah, mp_real tol, mp_real *totalerror)
//mp_real rk78()
/*
this routine performs one step of the integration procedure.
the initial condition (at,x) is changed by a new one corresponding
to the same orbit. the error is controlled by the threshold tol,
and an estimate of the error produced in the actual step is returned
as the value of the function.

parameters:
at:   time. input: time corresponding to the actual initial condition.
            output: new value corresponding to the new initial condition.
x:    position. same remarks as at.
ah:   time step (it can be modified by the routine according to the
      given threshold).
tol:  threshold to control the integration error.
hmin: minimun time step allowed.
hmax: maximum time step allowed.
n:    dimension of the system of odes.
deriv: function that returns the value of the vectorfield.

returned value: an estimate of the error produced in the actual step of
integration.
*/
{
   //if (n > neq) {printf("rk78: wrong dimension (%d and %d)\n",n,neq); exit(1);}

mp_real *alfa,*beta,*c8,*c7;

mp_real tpon,tol1,err,nor,kh,beth,h1;
mp_real hmax,hmin;
mp_real *x7,*x8,*xpon,*dx;
mp_real** k;
int n,i,j,m,l;
alfa=new mp_real[13];
beta=new mp_real[79];
c8=new mp_real[13];
c7=new mp_real[11];

x7=new mp_real[5];
x8=new mp_real[5];
xpon=new mp_real[5];
dx=new mp_real[5];
k=new mp_real*[13];
for (i=0;i<13;i++){
	k[i]=new mp_real[5];
}

obtain_constants(alfa,beta,c7,c8);
n=5;
hmin=mp_real("1e-20"); //Minimum step for the rk78
hmax=mp_real("1e-1"); //Maximum step for the rk78

//This is temporal:
/*
mp_real *x,*at,*xnext,params[5],delta,*ah,tol;
int rblock,lblock;
at=new mp_real;
xnext=new mp_real[5];
ah=new mp_real;
x=new mp_real[5];
params[0]=mp_real("2");
params[1]=mp_real("1");
params[2]=mp_real("1.5");
lblock=1;
rblock=1;
tol=mp_real("1e-5");
 */

do {
/*      
this is to compute the values of k
*/
      m=0;
      for (i=0; i<13; i++)
      {
         //tpon=*at+alfa[i]*(*ah);
         tpon=x[4]+alfa[i]*(*ah);

         for (j=0; j<n; j++ ) xpon[j]=x[j];
         for ( l=0; l<i; l++ )
         {
            ++m;
            beth=*ah*beta[m];
            for (j=0; j<n; j++) xpon[j] += beth*k[l][j];
         }
         //(*deriv)(tpon,xpon,n,dx);
         f(xpon,dx,params,delta,lblock,rblock);
         for (j=0; j<n; j++ ){
		//cout <<x[j]<<endl;
		k[i][j] = dx[j];
	}
      }
/*
      this is to compute the rk7 and rk8 predictions
*/
      err=nor=mp_real("0.e0");
      for (j=0; j<n; j++)
      {
         x7[j]=x8[j]=x[j];
         for (l=0; l<11; l++)
         {
            kh=*ah*k[l][j];
            x7[j] += kh*c7[l];
            x8[j] += kh*c8[l];
         }
         x8[j] += *ah*(c8[11]*k[11][j]+c8[12]*k[12][j]);
         err += abs(x8[j]-x7[j]);
         nor += abs(x8[j]);
      }
      err /= n;

/*
      next lines compute the new time step h
*/
      tol1=tol*(1+nor/mp_real("100"));
      if (err < tol1) err=MAX(err,tol1/mp_real("256"));
      h1=*ah;
      *ah*=mp_real("0.9")*pow(tol1/err,mp_real("0.125"));
      if (abs(*ah) < hmin ) *ah=hmin*SGN(*ah);
      if (abs(*ah) > hmax ) *ah=hmax*SGN(*ah);

   } while ((err >= tol1) && (abs (*ah) > hmin));
   //*at += h1;
   for (j=0; j<n; j++) xnext[j]=x8[j];

*totalerror=*totalerror+err;

//This is temporaL
/*
delete ah;
delete at;
delete [] xnext;
delete [] x;
*/
delete [] x7;
delete [] x8;
delete [] xpon;
delete [] dx;
for (i=0;i<13;i++){
delete [] k[i];
}
delete [] k;
delete [] alfa;
delete [] beta;
delete [] c8;
delete [] c7;
   //return (err);

}

mp_real SGN(mp_real x){

	if (x>=0){
		return mp_real("1");
	}
	else
		return -mp_real("1");
}

mp_real MAX(mp_real a,mp_real b){
	if (a<b){
		return b;
	}
	else
		return a;

}


void obtain_constants(mp_real *alfa,mp_real *beta,mp_real *c7,mp_real *c8){

alfa[0]=mp_real("0.e0");
alfa[1]=mp_real("2.e0")/mp_real("27.e0");
alfa[2]=mp_real("1.e0/9.e0");
alfa[3]=mp_real("1.e0")/mp_real("6.e0");
alfa[4]=mp_real("5.e0")/mp_real("12.e0");
alfa[5]=mp_real("0.5e0");
alfa[6]=mp_real("5.e0")/mp_real("6.e0");
alfa[7]=mp_real("1.e0")/mp_real("6.e0");
alfa[8]=mp_real("2.e0")/mp_real("3.e0");
alfa[9]=mp_real("1.e0")/mp_real("3.e0");
alfa[10]=mp_real("1.e0");
alfa[11]=mp_real("0.e0");
alfa[12]=mp_real("1.e0");

beta[0]=mp_real("0.e0");
beta[1]=mp_real("2.e0")/mp_real("27.e0");
beta[2]=mp_real("1.e0")/mp_real("36.e0");
beta[3]=mp_real("1.e0")/mp_real("12.e0");
beta[4]=mp_real("1.e0")/mp_real("24.e0");
beta[5]=mp_real("0.e0");
beta[6]=mp_real("1.e0")/mp_real("8.e0");
beta[7]=mp_real("5.e0")/mp_real("12.e0");
beta[8]=mp_real("0.e0");
beta[9]=mp_real("-25.e0")/mp_real("16.e0");
beta[10]=mp_real("25.e0")/mp_real("16.e0");
beta[11]=mp_real(".5e-1");
beta[12]=mp_real("0.e0");
beta[13]=mp_real("0.e0");
beta[14]=mp_real(".25e0");
beta[15]=mp_real(".2e0");
beta[16]=mp_real("-25.e0")/mp_real("108.e0");
beta[17]=mp_real("0.e0");
beta[18]=mp_real("0.e0");
beta[19]=mp_real("125.e0")/mp_real("108.e0");
beta[20]=mp_real("-65.e0")/mp_real("27.e0");
beta[21]=mp_real("125.e0")/mp_real("54.e0");
beta[22]=mp_real("31.e0")/mp_real("300.e0");
beta[23]=mp_real("0.e0");
beta[24]=mp_real("0.e0");
beta[25]=mp_real("0.e0");
beta[26]=mp_real("61.e0")/mp_real("225.e0");
beta[27]=mp_real("-2.e0")/mp_real("9.e0");
beta[28]=mp_real("13.e0")/mp_real("900.e0");
beta[29]=mp_real("2.e0");
beta[30]=mp_real("0.e0");
beta[31]=mp_real("0.e0");
beta[32]=mp_real("-53.e0")/mp_real("6.e0");
beta[33]=mp_real("704.e0")/mp_real("45.e0");
beta[34]=mp_real("-107.e0")/mp_real("9.e0");
beta[35]=mp_real("67.e0")/mp_real("90.e0");
beta[36]=mp_real("3.e0");
beta[37]=mp_real("-91.e0")/mp_real("108.e0");
beta[38]=mp_real("0.e0");
beta[39]=mp_real("0.e0");
beta[40]=mp_real("23.e0")/mp_real("108.e0");
beta[41]=mp_real("-976.e0")/mp_real("135.e0");
beta[42]=mp_real("311.e0")/mp_real("54.e0");
beta[43]=mp_real("-19.e0")/mp_real("60.e0");
beta[44]=mp_real("17.e0")/mp_real("6.e0");
beta[45]=mp_real("-1.e0")/mp_real("12.e0");
beta[46]=mp_real("2383.e0")/mp_real("4100.e0");
beta[47]=mp_real("0.e0");
beta[48]=mp_real("0.e0");
beta[49]=mp_real("-341.e0")/mp_real("164.e0");
beta[50]=mp_real("4496.e0")/mp_real("1025.e0");
beta[51]=mp_real("-301.e0")/mp_real("82.e0");
beta[52]=mp_real("2133.e0")/mp_real("4100.e0");
beta[53]=mp_real("45.e0")/mp_real("82.e0");
beta[54]=mp_real("45.e0")/mp_real("164.e0");
beta[55]=mp_real("18.e0")/mp_real("41.e0");
beta[56]=mp_real("3.e0")/mp_real("205.e0");
beta[57]=mp_real("0.e0");
beta[58]=mp_real("0.e0");
beta[59]=mp_real("0.e0");
beta[60]=mp_real("0.e0");
beta[61]=mp_real("-6.e0")/mp_real("41.e0");
beta[62]=mp_real("-3.e0")/mp_real("205.e0");
beta[63]=mp_real("-3.e0")/mp_real("41.e0");
beta[64]=mp_real("3.e0")/mp_real("41.e0");
beta[65]=mp_real("6.e0")/mp_real("41.e0");
beta[66]=mp_real("0.e0");
beta[67]=mp_real("-1777.e0")/mp_real("4100.e0");
beta[68]=mp_real("0.e0");
beta[69]=mp_real("0.e0");
beta[70]=mp_real("-341.e0")/mp_real("164.e0");
beta[71]=mp_real("4496.e0")/mp_real("1025.e0");
beta[72]=mp_real("-289.e0")/mp_real("82.e0");
beta[73]=mp_real("2193.e0")/mp_real("4100.e0");
beta[74]=mp_real("51.e0")/mp_real("82.e0");
beta[75]=mp_real("33.e0")/mp_real("164.e0");
beta[76]=mp_real("12.e0")/mp_real("41.e0");
beta[77]=mp_real("0.e0");
beta[78]=mp_real("1.e0");

c7[0]=mp_real("41.e0")/mp_real("840.e0");
c7[1]=mp_real("0.e0");
c7[2]=mp_real("0.e0");
c7[3]=mp_real("0.e0");
c7[4]=mp_real("0.e0");
c7[5]=mp_real("34.e0")/105.e0;
c7[6]=mp_real("9.e0")/35.e0;
c7[7]=mp_real("9.e0")/mp_real("35.e0");
c7[8]=mp_real("9.e0")/mp_real("280.e0");
c7[9]=mp_real("9.e0")/mp_real("280.e0");
c7[10]=mp_real("41.e0")/mp_real("840.e0");

c8[0]=mp_real("0.e0");
c8[1]=mp_real("0.e0");
c8[2]=mp_real("0.e0");
c8[3]=mp_real("0.e0");
c8[4]=mp_real("0.e0");
c8[5]=mp_real("34.e0")/mp_real("105.e0");
c8[6]=mp_real("9.e0")/mp_real("35.e0");
c8[7]=mp_real("9.e0")/mp_real("35.e0");
c8[8]=mp_real("9.e0")/mp_real("280.e0");
c8[9]=mp_real("9.e0")/mp_real("280.e0");
c8[10]=mp_real("0.e0");
c8[11]=mp_real("41.e0")/mp_real("840.e0");
c8[12]=mp_real("41.e0")/mp_real("840.e0");

}
