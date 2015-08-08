#ifndef _CORE_MC_FUNC_H_
#define _CORE_MC_FUNC_H_


#define maxClusters 10000
typedef struct {
	double x;
	double y;
	double z;
	int cluster;
}particle;//refers to an individual particle

typedef struct{
	int numClusters;
	particle clusters[maxClusters];
}clusters;//refers to a set of clusters with associated centres of masses

void Move_Particles(particle *particles,double *randoms,long int *currentRandom, int *accepted, int *rejected,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist);//moves each of the particles ones
void Init_Randoms(double *randoms,long int *currentRandom);
void FetchRands(int nrand,double *randoms);
double PRand(double *randoms,long int *currentRandom);//returns a random value between 0 and 1
int Accept(double BF_old,double BF_new,double *randoms,long int *currentRandom);//decides if a change is accepted based on boltzmann factors(1 accepted, 0 rejected)
double BFactor(double E_r,double temp);//calculates Boltzmann factor given potenial
double Potential(particle x,particle *particles,double xMax,double yMax,double zMax,double truncation,int numParticles);//calculates potential at a given point
double PairPotential(double r);//gives the potential from distance r
double Distance(particle a,particle b);//gives the distance between two points
double SDistance(particle a,particle b,double xMax,double yMax,double zMax);//gives the shortest distance between two points assuming periodic boundary conditions
double mod(double a,double b);//returns a%b in a way that works better for us
#endif