#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include"coreMCfunc.h"

void Move_Particles(particle *particles,double *randoms,long int *currentRandom,int *accepted, int *rejected,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist){
	int j;
	particle move;
	double bf1,bf2;
	for(j=0;j<numParticles;j++){
		//set trail move
		move.x = (particles+j)->x+2*(PRand(randoms,currentRandom)-0.5)*moveDist;
		move.y = (particles+j)->y+2*(PRand(randoms,currentRandom)-0.5)*moveDist;
		move.z = (particles+j)->z+2*(PRand(randoms,currentRandom)-0.5)*moveDist;

		//implement boundary conditions
		move.x = mod(move.x,xMax);
		move.y = mod(move.y,yMax);
		move.z = mod(move.z,zMax);
		//compare boltzman factors and decide if move is accepted
		bf1 = BFactor(Potential(*(particles+j),particles,xMax,yMax,zMax,truncation,numParticles) - PairPotential(0),temp);//second term accounts for the fact that a particle will not interact with itself
		bf2 = BFactor(Potential(move,particles,xMax,yMax,zMax,truncation,numParticles) - PairPotential(SDistance(*(particles+j),move,xMax,yMax,zMax)),temp);//second term removes contribution from particle you are moving 
		if(Accept(bf1,bf2,randoms,currentRandom)){
			*(particles+j) = move;
			*(accepted) += 1;
		}
		else{
			*(rejected) += 1;
		}
	}
	return;
}

void Init_Randoms(double *randoms,long int *currentRandom){
	int seed,nrand = 100000;
	*(currentRandom) = 0;
	# ifdef DEBUG
  	seed = 1;
	# else
	time((time_t *)&seed);
	# endif
	srandom(seed);
	FetchRands(nrand,randoms);
	return;
}

void FetchRands(int nrand,double *randoms){
	int i;
	for(i=0;i<nrand;i++){
		*(randoms+i) = (double)random()/(double)RAND_MAX;
	}
}

double PRand(double *randoms,long int *currentRandom){
	if(*(currentRandom) >= 100000){
		FetchRands(100000,randoms);
		*(currentRandom) = 0;
	}
	double out = randoms[*(currentRandom)];
	*(currentRandom) += 1;
	return out;
}

int Accept(double BF_old,double BF_new,double *randoms,long int *currentRandom){
	if(PRand(randoms,currentRandom)<(BF_new/BF_old)){
		return 1;
	}
	else{
		return 0;
	}
}

double BFactor(double E_r,double temp){
	double ret = exp(-E_r/(temp));
	return ret;
}

double Potential(particle x,particle *particles,double xMax,double yMax,double zMax,double truncation,int numParticles){
	int i;
	double d;//distance between particles(temporary)
	double out = 0;//running total of potential
	double truncationDist = truncation;
	for(i=0;i<numParticles;i++){
		//finds shortest distance to an instance of particles[i]
		d = SDistance(*(particles+i),x,xMax,yMax,zMax);
		//adds potential to sum if d is bellow trunctation cutoff
		if(d<truncationDist){
			out += PairPotential(d);
		}
	}
	return out;
}




double PairPotential(double r){
	double ret =  exp(-pow((r),4));
	return ret;
}

double Distance(particle a,particle b){
	double Delta_x = a.x-b.x;
	double Delta_y = a.y-b.y;
	double Delta_z = a.z-b.z;
	return sqrt((Delta_x*Delta_x)+(Delta_y*Delta_y)+(Delta_z*Delta_z));
}

double SDistance(particle a,particle b,double xMax,double yMax,double zMax){
	double dx,dy,dz,d;
	dx = fabs(a.x-b.x)-(floor((fabs(a.x-b.x)/xMax)+0.5)*xMax);
	dy = fabs(a.y-b.y)-(floor((fabs(a.y-b.y)/yMax)+0.5)*yMax);
	dz = fabs(a.z-b.z)-(floor((fabs(a.z-b.z)/zMax)+0.5)*zMax);
	d = sqrt((dx*dx)+(dy*dy)+(dz*dz));
	return d;
}

double mod(double a,double b){//cant deal with negative b
	double out;
	if(a>=0){
		out = (a-(floor(a/b)*b));
	}
	else{
		out = (a-(-floor(a/-b)*b));//to deal with floor working differently on negative numbers
	}
	if(out<0){
		return out+b;
	}
	return out;
}