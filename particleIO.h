#ifndef _PARTICLE_IO_H_
#define _PARTICLE_IO_H_
void Populate_Random(particle *particles,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles);//adds numParticles particles randomly to the array particles
void Populate_Cubic(particle *particles,double xMax,int numParticles);//adds numparticles to the array particles in a simple cubic stucture (use cube values of n)
void Output_print(particle *particles,int numToOutput);//prints the current state of the system
void Populate_Cluster_fcc(particle *particles,double a,double spread,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles);//add numparticles to particles in a fcc structure of lattice constant a, adding multiple occupancy after filling once
void Populate_Cluster_bcc(particle *particles,double a,double spread,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles);//add numparticles to particles in a bcc structure of lattice constant a, adding multiple occupancy after filling once
void Output_jmol(particle *particles,char *name,double scaling,int numToOutput,double xMax,double yMax,double zMax);//outputs the current state of the system as a jmol file with filename name.xyz(scaling factor is used to make the file easier to visualise)
void Input_jmol(particle *particles,char *name,double scaling,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles);//reads from name.xyz to set the current state to that of a past simulation(must know scaling factor used when output)
#endif