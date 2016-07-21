#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int main (int argc, char *argv[])
{

double tend1,tend2,tini1,tini2;
double avf,avb;
double x,y,u,v,s,s0,prevs;
double U,aux;
FILE *forwards,*backwards,*fout1,*fout2;


forwards=fopen("heteroclinic_forwards.dat","r");
backwards=fopen("heteroclinic_backwards.dat","r");
fout1=fopen("averaged_forwards.dat","w");
fout2=fopen("averaged_backwards.dat","w");

//If tini is 0, we integrate from s=s0. If it is not, we start the windows at
//s=tini.
//Careful, the first value of s is not s0 but s0+zeta^* because we integrate the
//flow with initial condition z_0^* and not z^*.

tini1=strtod(argv[1],NULL);
tend1=strtod(argv[2],NULL);
tini2=strtod(argv[3],NULL);
tend2=strtod(argv[4],NULL);

avf=0;
avb=0;
aux=0;


fscanf(forwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s0);
s=s0;
if (tini1!=0){
	while (s<tini1){
		fscanf(forwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s);
	}
}

while (s<tend1 && !feof(forwards)){
	prevs=s;
	fscanf(forwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s);
	U=pow(v,2)/2-pow(u,2)/2+fabs(u);
	aux=aux+U*(s-prevs);
	avf=aux/(s-s0);
	fprintf(fout1,"%lf %lf\n",s-s0,avf);
}


aux=0;

fscanf(backwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s0);
s=s0;
if (tini2!=0){
	while (s>tini2){
		fscanf(backwards,"%lf %lf %lf %lf %lf",&x,&y,&u,&v,&s);
	}
}
//we assume that tend2 is negative!
while (s>tend2 && !feof(backwards)){
	prevs=s;
	fscanf(backwards,"%lf %lf %lf %lf %lf\n",&x,&y,&u,&v,&s);
	U=pow(v,2)/2-pow(u,2)/2+fabs(u);
	aux=aux+U*(prevs-s);
	avb=aux/(s0-s);
	fprintf(fout2,"%lf %lf\n",s0-s,avb);
}

fclose(forwards);
fclose(backwards);
fclose(fout1);
fclose(fout2);

}
