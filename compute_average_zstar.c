#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int main (int argc, char *argv[])
{

double tend1,tend2;
double avf,avb;
double x,y,u,v,s,s0,prevs;
double U,auxf,auxb;
FILE *forwards,*backwards,*fout1,*fout2;


forwards=fopen("heteroclinic_forwards.dat","r");
backwards=fopen("heteroclinic_backwards.dat","r");
fout1=fopen("averaged_forwards.dat","w");
fout2=fopen("averaged_backwards.dat","w");

//If tini is 0, we integrate from s=s0. If it is not, we start the windows at
//s=tini.
//Careful, the first value of s is not s0 but s0+zeta^* because we integrate the
//flow with initial condition z_0^* and not z^*.
//In this version we compute the average integrating from z^*. Assuming s0=0,
//z^* corresponds to s=0, which assume is located in heteroclinic_backwards.dat,
//because we are dealing with positive zeros of the Melnikov function.

tend1=strtod(argv[1],NULL);
tend2=strtod(argv[2],NULL);

avf=0;
avb=0;
auxf=0;
auxb=0;


//We start with heteroclinic_backwards.dat
fscanf(backwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s0);
s=s0;

while (s>0 && !feof(backwards)){
	//While s>0 this goes for the forwards integration.
	prevs=s;
	fscanf(backwards,"%lf %lf %lf %lf %lf\n",&x,&y,&u,&v,&s);
	U=pow(v,2)/2-pow(u,2)/2+fabs(u);
	auxf=auxf+U*(prevs-s);
	avf=auxf/(s0-s);
	fprintf(fout1,"%lf %lf\n",s0-s,avf);
}

//Now s<0 and this goes to backwards integration.
//we assume that tend2 is negative!
while (s>tend2 && !feof(backwards)){
	prevs=s;
	fscanf(backwards,"%lf %lf %lf %lf %lf\n",&x,&y,&u,&v,&s);
	U=pow(v,2)/2-pow(u,2)/2+fabs(u);
	auxb=auxb+U*(prevs-s);
	avb=auxb/(-s);
	fprintf(fout2,"%lf %lf\n",-s,avb);
}



fscanf(forwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s);

while (s<tend1 && !feof(forwards)){
	prevs=s;
	fscanf(forwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s);
	U=pow(v,2)/2-pow(u,2)/2+fabs(u);
	auxf=auxf+U*(s-prevs);
	avf=auxf/(s);
	fprintf(fout1,"%lf %lf\n",s,avf);
}



fclose(forwards);
fclose(backwards);
fclose(fout1);
fclose(fout2);

}
