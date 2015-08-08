#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_
void Output_gr(char *name,int numToOutput,particle *particles,double xMax,double yMax,double zMax);//outputs g(r) to name.graph
clusters ID_Clusters(particle *particles,double *avgclusterSize,int includingSmallClusters,double xMax,double yMax,double zMax,int numParticles,double clusterDist);//identifies the particles into clusters and prints the average cluster size as well as setting the variable pointed to
void Output_Cluster_Distribution(char *name,clusters clusters,particle *particles,double avgClusterSize,double xMax,double yMax,double zMax,int numParticles);//outputs the distribution of cluster sizes to name.dist
void Output_100avg_AvgCluster(particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,double clusterDist);//takes the average cluster size of the system over 100 MC cycles & give the mean and std deviation
void Output_100avg_ClusterDist(char *name,particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,double clusterDist);//takes cluster dist over 100 cycles, then outputs avg, min &max estimates based on varience
void Output_100avg_gr(char *name,particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected, int numToOutput,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist);//takes g(r) over 100 cycles and then outputs avg, min & max based on varience.
void Output_100avg_Cluster_gr(char *name,particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected, int numToOutput,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,double clusterDist);//takes cluster g(r) over 100 cycles and then outputs avg, min & max based on varience.
#endif