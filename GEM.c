/****************************************************************************************
*program to simulate a generalised exponential potential model				*
*											*
*this is a (slightly tidied up version of) the program written by       		*
*myself (stu.wood.94@googlemail.com) and Jack Porter as part of our			*
*final year project.									*
*											*
*header files also contain descriptions of the purpose of each function contained within*
****************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include"coreMCfunc.h"//contains functions essential to the running monte carlo cycles on the system
#include"particleIO.h"//contains functions linked to input and output and population of particle states
#include"analysis.h"//contains functions linked to analyising a configuration of particles

//define DEBUG //include this line to give same random numbers each time

void mainCycle(particle * particles,double * randoms,long int * currentRandom,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,int cycles);//this is the meat of the simulation, it runs the MC cycles to move the system into a low energy state.

int main(int argc,char*argv[]){
  
	//parameters
	int numParticles = 450*8;//the number of particles in the first simulation, density is this number/100
	int repeats = 1;
	int particleIncrement = 200;//the amount that the number of particles increases by each loop
	int particleLoops = 1;//the number of different densities examined
	double temp = 4;//temperaure of first simulation
	double tempIncrement = 0.2;//the amount that the temp increases by each loop
	int tempLoops = 1;//the number of different temps examined
	double xMax = 2*4.641588834;//   ^
	double yMax = 2*4.641588834;//   |  theses are the sizes for the boundaries
	double zMax = 2*4.641588834;//   v
	double moveDist = 0.4;//maximum translational move
	int cycles =10;//number of moves of each particle(set at >= 10 if running any cycles
	double truncation = 2;//the distance at whcih potential is ignored-measured in multiples of sigma
	double clusterDist = 0.25;//sets the maximum distance at which particles are considered to be in the same cluster
	
	//initializing random number generator
	double *randoms = malloc(sizeof(double)*100000);//used to store random numbers & avoid repeatedly calling ranvec
	long int currentRandom = 0;//indicates how far through randoms[] we are
	Init_Randoms(randoms,&currentRandom);

	//declare variables and allocate memory
	int k,j,l;
	particle *particles = malloc(sizeof(particle)*(numParticles+particleLoops*particleIncrement));
	clusters *bigClusters = malloc(sizeof(clusters));
	clusters *allClusters = malloc(sizeof(clusters));
	double startTemp = temp;
	double startParticles = numParticles;
	double avgClusterSize;

	//start simulations
	//these loops allow for multiple simulations in different conditions with one execution
	for(l=0;l<repeats;l++){
		numParticles = startParticles;
		for(k=0;k<particleLoops;k++){
			temp = startTemp;
			for(j=0;j<tempLoops;j++){

				//poulate the system
				Populate_Random(particles,randoms,&currentRandom,xMax,yMax,zMax,numParticles);

				//main cycle
				mainCycle(particles,randoms,&currentRandom,xMax,yMax,zMax,truncation,numParticles,temp,moveDist,cycles);
				//output measures for final state*/
				*(bigClusters) = ID_Clusters(particles,&avgClusterSize,0,xMax,yMax,zMax,numParticles,clusterDist);
				*(allClusters) = ID_Clusters(particles,&avgClusterSize,1,xMax,yMax,zMax,numParticles,clusterDist);
				Output_Cluster_Distribution("Clusters",*(allClusters),particles,avgClusterSize,xMax,yMax,zMax,numParticles);
				Output_jmol(bigClusters->clusters,"Clusters",5,bigClusters->numClusters,xMax,yMax,zMax);
				Output_gr("Clusters",bigClusters->numClusters,bigClusters->clusters,xMax,yMax,zMax);
				Output_jmol(particles,"Particles",5,numParticles,xMax,yMax,zMax);
				Output_gr("Particles",numParticles,particles,xMax,yMax,zMax);
			
				//update parameters for next simulation
				temp += tempIncrement;
			}
			numParticles += particleIncrement;
		}
	}
	
	//cleanup
	free(bigClusters);
	free(allClusters);
	free(particles);
	free(randoms);
	return 0;
}

void mainCycle(particle * particles,double * randoms,long int * currentRandom,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,int cycles){
				int i;
				int accepted = 0;//counts accepted moves
				int rejected = 0;//counts rejected moves
				if(xMax<(2*truncation)){//checks that box is at least double truncation distance
					printf("warning:box is too small for current parameters\n");
				}
				for(i=0;i<cycles;i++){
					Move_Particles(particles,randoms,currentRandom,&accepted,&rejected,xMax,yMax,zMax,truncation,numParticles,temp,moveDist);
					//prints progress
					if(i%(cycles/10)==0){
						printf("d%.2lfT%.2lfm%.2lfc%d: %d%%   acceptance rate:%lf\n",(double)numParticles/(xMax*yMax*zMax),temp,moveDist,cycles,(i*100)/cycles,(double)accepted/(accepted+rejected));
					}
				}
				printf("acceptance rate: %lf\n",(double)accepted/(double)(accepted+rejected));
				return;
}