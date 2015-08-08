#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include"coreMCfunc.h"
#include"analysis.h"

void Output_gr(char *name,int numToOutput,particle *particles,double xMax,double yMax,double zMax){
	int xMaxi = ceil(xMax);
	int i,j;
	int dr = 50;//the resiprical of the spacing
	//initialises the count array
	long int num_at_dr[xMaxi*dr];
	for(i=0;i<xMaxi*dr;i++){
		num_at_dr[i] = 0;
	}
	int distances[numToOutput];
	for(j=0;j<numToOutput;j++){
		//store the closest distances of each of the particles
		for(i=0;i<numToOutput;i++){
			distances[i] = floor(dr*SDistance(*(particles+i),*(particles+j),xMax,yMax,zMax));//this line discretises the distances into multiples of 1/dr after finding them, leaving them as distance*dr
		}
		//count the number of particles in each distance value.
		for(i=0;i<numToOutput;i++){
			num_at_dr[distances[i]]++;
		}
	}
	//find corresponding g(r) for each distance value up to xMax
	double gr[xMaxi*dr];
	gr[0] = 0;
	for(i=1;i<xMaxi*dr;i++){
		gr[i] = ((double)num_at_dr[i]*xMax*yMax*zMax*(double)dr*(double)(pow(dr,2)))/((double)4*(double)(i)*(double)(i)*(double)numToOutput*(double)numToOutput*M_PI);//the dr^2 is to make normalising work, not sure why we need it
	}
	FILE *out;
	char lname[strlen(name)+6];
	strcpy(lname,name);
	strcat(lname,".graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	char nameEnd[9];
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,gr[i]);
	}
	fclose(out);
	printf("g(r) saved to %s\n",lname);
	return;
}

clusters ID_Clusters(particle *particles,double *avg,int includingSmallClusters,double xMax,double yMax,double zMax,int numParticles,double clusterDist){
	int i,j;
	int numClusters = 0;
	int stack[numParticles];
	int numInStack = 0;
	//first we reset any previous IDing of clusters
	for(i=0;i<numParticles;i++){
		(particles+i)->cluster = 0;
	}
	//now we go about retagging them
	for(i=0;i<numParticles;i++){
		if((particles+i)->cluster==0){
			//found a new cluster
			numClusters++;
			if(numClusters>maxClusters){
				printf("too many clusters, aborting ID\n");
				clusters ret;
				return ret;
			}
			(particles+i)->cluster = numClusters;//give new cluster an ID
			stack[0] = i;//add the first element of the cluster to the stack
			numInStack = 1;
			//look for others in cluster
			while(numInStack>0){
				for(j=0;j<numParticles;j++){
					if(j!=stack[0]&&SDistance(*(particles+j),*(particles+stack[0]),xMax,yMax,zMax)<clusterDist&&(particles+j)->cluster==0){
						//new member of cluster found
						(particles+j)->cluster = numClusters;//set it's ID
						stack[numInStack] = j;//add it to the stack
						numInStack++;
					}
				}
				//finished looking around this particle, take it off the stack
				stack[0] = stack[numInStack-1];
				numInStack--;
			}
		}
	}
	//now we find the aveage of each cluster
	particle totals[numClusters];
	int sizeOfClusters[numClusters];
	clusters allClusters;
	clusters bigClusters;
	//initialising
	for(i=0;i<numClusters;i++){
		sizeOfClusters[i] = 0;
		totals[i].x = 0;
		totals[i].y = 0;
		totals[i].z = 0;
	}
	//Total the components for each particle by cluster
	for(i=0;i<numParticles;i++){
		sizeOfClusters[((particles+i)->cluster)-1]++;
		totals[((particles+i)->cluster)-1].x += (particles+i)->x;
		totals[((particles+i)->cluster)-1].y += (particles+i)->y;
		totals[((particles+i)->cluster)-1].z += (particles+i)->z;
	}
	//use totals to find cluster centre points
	for(i=0;i<numClusters;i++){
		allClusters.clusters[i].x = (totals[i].x/(double)sizeOfClusters[i]);
		allClusters.clusters[i].y = (totals[i].y/(double)sizeOfClusters[i]);
		allClusters.clusters[i].z = (totals[i].z/(double)sizeOfClusters[i]);
	}
	allClusters.numClusters = numClusters;
	double avgClusterSize = (double)numParticles/numClusters;
	*avg=avgClusterSize;
	printf("Average cluster size: %lf\n",avgClusterSize);
	if(includingSmallClusters){
		return allClusters;
	}
	//now we copy over to bigClusters, removing any clustser of size<half avg
	bigClusters.numClusters = allClusters.numClusters;
	j = 0;
	for(i=0;i<numClusters;i++){
		if(sizeOfClusters[i]<(avgClusterSize/2)){
			bigClusters.numClusters -= 1;
			continue;
		}
		else{
			bigClusters.clusters[j] = allClusters.clusters[i];
			j++;
		}
	}
	return bigClusters;
}

void Output_Cluster_Distribution(char *name,clusters clusters,particle *particles,double avgClusterSize,double xMax,double yMax,double zMax,int numParticles){
	int maxBin = ceil(30);
	int i,j;
	int num_of_size[maxBin];
	int size_of_cluster[clusters.numClusters];
	//initialise the counting array
	for(i=0;i<maxBin;i++){
		num_of_size[i] = 0;
	}
	//find the size of each cluster
	for(i=1;i<clusters.numClusters;i++){
		size_of_cluster[i] = 0;
		for(j=0;j<numParticles;j++){
			if(((particles+j)->cluster) == i){
				size_of_cluster[i]++;
			}
		}
	}
	//counts the number of clusters of each size
	for(i=0;i<clusters.numClusters;i++){
		num_of_size[size_of_cluster[i]]++;
	}
	//outputs our data
	FILE *out;
	char lname[strlen(name)+6];
	strcpy(lname,name);
	strcat(lname,".dist");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	char nameEnd[9];
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"%d.dist",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<maxBin;i++){
		fprintf(out,"%d	%d\n",i,num_of_size[i]);
	}
	fclose(out);
	printf("cluster distribution saved to %s\n",lname);
	return;
}

void Output_100avg_AvgCluster(particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,double clusterDist){
	double reads[100];
	int i;
	//first gather our data set
	for(i=0;i<100;i++){
		ID_Clusters(particles,(reads+i),1,xMax,yMax,zMax,numParticles,clusterDist);
		Move_Particles(particles, randoms, currentRandom, accepted, rejected,xMax,yMax,zMax,truncation,numParticles,temp,moveDist);
		printf("%d   ",i);
	}
	//now find the mean
	double average = 0;
	double stdDev;
	for(i=0;i<100;i++){
		average += reads[i];
	}
	average = average/100;
	//now find the standard deviation
	for(i=0;i<100;i++){
		stdDev += (reads[i]-average)*(reads[i]-average);
	}
	stdDev = sqrt(stdDev/100);
	printf("mean = %lf\nstdDev = %lf\n", average, stdDev);
	return;
}

void Output_100avg_ClusterDist(char *name,particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,double clusterDist){
	int reads_num_of_size[30][100];
	clusters clusters;
	int num_of_size[30];
	int *size_of_cluster = malloc(sizeof(int)*maxClusters);
	int i,j,k;
	double dummy;
//gather data
	for(k=0;k<100;k++){
		clusters = ID_Clusters(particles,&dummy, 1,xMax,yMax,zMax,numParticles,clusterDist);
		//initialise the counting array
		for(i=0;i<30;i++){
			num_of_size[i] = 0;
		}
		//find the size of each cluster(unneccesarily slow)
		for(i=1;i<clusters.numClusters;i++){
			size_of_cluster[i] = 0;
			for(j=0;j<numParticles;j++){
				if(((particles+j)->cluster) == i){
					size_of_cluster[i]++;
				}
			}
		}
		//counts the number of clusters of each size
		for(i=0;i<clusters.numClusters;i++){
			num_of_size[size_of_cluster[i]]++;
		}
		//writes values onto reads
		for(i=0;i<30;i++){
			reads_num_of_size[i][k] = num_of_size[i];
		}
		//move the system on
		Move_Particles(particles, randoms, currentRandom, accepted, rejected,xMax,yMax,zMax,truncation,numParticles,temp,moveDist);
		printf("%d   ",k);
	}
//take the averages
	double average[30];
	double stdDev[30];
	for(i=0;i<30;i++){
		average[i] = 0;
		for(j=0;j<100;j++){
			average[i] += reads_num_of_size[i][j];
		}
		average[i] = average[i] / 100;
		stdDev[i] = 0;
		for(j=0;j<100;j++){
			stdDev[i] += (reads_num_of_size[i][j]-average[i])*(reads_num_of_size[i][j]-average[i]);
		}
		stdDev[i] = sqrt(stdDev[i]/100);
	}
	//now we find the max and min possible values of the distribution based on stdDev
	double maxEst[30];
	double minEst[30];
	for(i=0;i<30;i++){
		maxEst[i] = average[i] + (stdDev[i]/10);
		minEst[i] = average[i] + (stdDev[i]/10);
	}
//now we must output our estimates
//average
	FILE *out;
	char lname[strlen(name)+9];
	strcpy(lname,name);
	strcat(lname,"est.dist");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	char nameEnd[12];
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"est%d.dist",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<30;i++){
		fprintf(out,"%d	%lf\n",i,average[i]);
	}
	fclose(out);
	printf("average cluster distribution saved to %s\n",lname);
//max
	strcpy(lname,name);
	strcat(lname,"max.dist");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"max%d.dist",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<30;i++){
		fprintf(out,"%d	%lf\n",i,maxEst[i]);
	}
	fclose(out);
	printf("maximum estimated cluster distribution saved to %s\n",lname);
//min
	strcpy(lname,name);
	strcat(lname,"min.dist");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"min%d.dist",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<30;i++){
		fprintf(out,"%d	%lf\n",i,minEst[i]);
	}
	fclose(out);
	printf("maximum estimated cluster distribution saved to %s\n",lname);
	free(size_of_cluster);
	return;
}

void Output_100avg_gr(char *name, particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected, int numToOutput,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist){
	int xMaxi = ceil(xMax);
	int i,j,k;
	int dr = 50;//the resiprical of the spacing
	int distances[numToOutput];
	long int num_at_dr[xMaxi*dr];
	double read_gr[xMaxi*dr][100];
//first we need our sample
	for(k=0;k<100;k++){
		//initialises the count array
		for(i=0;i<xMaxi*dr;i++){
			num_at_dr[i] = 0;
		}
		for(j=0;j<numToOutput;j++){
			//store the closest distances of each of the particles
			for(i=0;i<numToOutput;i++){	
				distances[i] = floor(dr*SDistance(*(particles+i),*(particles+j),xMax,yMax,zMax));//this line discretises the distances into multiples of 1/dr after finding them, leaving them as distance*dr
			}
			//count the number of particles in each distance value.
			for(i=0;i<numToOutput;i++){
				num_at_dr[distances[i]]++;
			}
		}
		//find corresponding g(r) for each distance value up to xMax
		read_gr[0][k] = 0;
		for(i=1;i<xMaxi*dr;i++){
			read_gr[i][k] = ((double)num_at_dr[i]*xMax*yMax*zMax*(double)dr*(double)(pow(dr,2)))/((double)4*(double)(i)*(double)(i)*(double)numToOutput*(double)numToOutput*M_PI);//the dr^2 is to make normalising work, not sure why we need it
		}
		Move_Particles(particles,randoms,currentRandom,accepted,rejected,xMax,yMax,zMax,truncation,numParticles,temp,moveDist);
		printf("%d   ",k);
	}
//now we take averages
	double average[xMaxi*dr];
	double stdDev[xMaxi*dr];
	for(i=0;i<xMaxi*dr;i++){
		average[i] = 0;
		for(j=0;j<100;j++){
			average[i] += read_gr[i][j];
		}
		average[i] = average[i]/100;
		stdDev[i] = 0;
		for(j=0;j<100;j++){
			stdDev[i] += (read_gr[i][j]-average[i])*(read_gr[i][j]-average[i]);
		}
		stdDev[i] = sqrt(stdDev[i]/100);
	}
	//now we get min and max estimates
	double minEst[xMaxi*dr];
	double maxEst[xMaxi*dr];
	for(i=0;i<xMaxi*dr;i++){
		minEst[i] = average[i] - (stdDev[i]/10);
		maxEst[i] = average[i] + (stdDev[i]/10);
	}
//time to output
//average
	FILE *out;
	char lname[strlen(name)+9];
	strcpy(lname,name);
	strcat(lname,"avg.graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	char nameEnd[12];
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"avg%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,average[i]);
	}
	fclose(out);
	printf("average g(r) saved to %s\n",lname);
//min
	strcpy(lname,name);
	strcat(lname,"min.graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"min%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,minEst[i]);
	}
	fclose(out);
	printf("estimated minimum g(r) saved to %s\n",lname);
//max
	strcpy(lname,name);
	strcat(lname,"max.graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"max%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,maxEst[i]);
	}
	fclose(out);
	printf("estimated maximum g(r) saved to %s\n",lname);
	return;
}

void Output_100avg_Cluster_gr(char *name, particle *particles, double *randoms, long int *currentRandom,int *accepted, int *rejected, int numToOutput,double xMax,double yMax,double zMax,double truncation,int numParticles,double temp,double moveDist,double clusterDist){
	double dummy;
	clusters clusters;
	int xMaxi = ceil(xMax);
	int i,j,k;
	int dr = 50;//the resiprical of the spacing
	int *distances = malloc(sizeof(int)*numParticles);
	long int num_at_dr[xMaxi*dr];
	double read_gr[xMaxi*dr][100];
//first we need our sample
	for(k=0;k<100;k++){
		clusters = ID_Clusters(particles,&dummy,1,xMax,yMax,zMax,numParticles,clusterDist);
		//initialises the count array
		for(i=0;i<xMaxi*dr;i++){
			num_at_dr[i] = 0;
		}
		for(j=0;j<clusters.numClusters;j++){
			//store the closest distances of each of the particles
			for(i=0;i<clusters.numClusters;i++){	
				distances[i] = floor(dr*SDistance(*(clusters.clusters+i),*(clusters.clusters+j),xMax,yMax,zMax));//this line discretises the distances into multiples of 1/dr after finding them, leaving them as distance*dr
			}
			//count the number of particles in each distance value.
			for(i=0;i<clusters.numClusters;i++){
				num_at_dr[distances[i]]++;
			}
		}
		//find corresponding g(r) for each distance value up to xMax
		read_gr[0][k] = 0;
		for(i=1;i<xMaxi*dr;i++){
			read_gr[i][k] = ((double)num_at_dr[i]*xMax*yMax*zMax*(double)dr*(double)(pow(dr,2)))/((double)4*(double)(i)*(double)(i)*(double)clusters.numClusters*(double)clusters.numClusters*M_PI);//the dr^2 is to make normalising work, not sure why we need it
		}
		Move_Particles(particles,randoms,currentRandom,accepted,rejected,xMax,yMax,zMax,truncation,numParticles,temp,moveDist);
		printf("%d   ",k);
	}
//now we take averages
	double average[xMaxi*dr];
	double stdDev[xMaxi*dr];
	for(i=0;i<xMaxi*dr;i++){
		average[i] = 0;
		for(j=0;j<100;j++){
			average[i] += read_gr[i][j];
		}
		average[i] = average[i]/100;
		stdDev[i] = 0;
		for(j=0;j<100;j++){
			stdDev[i] += (read_gr[i][j]-average[i])*(read_gr[i][j]-average[i]);
		}
		stdDev[i] = sqrt(stdDev[i]/100);
	}
	//now we get min and max estimates
	double minEst[xMaxi*dr];
	double maxEst[xMaxi*dr];
	for(i=0;i<xMaxi*dr;i++){
		minEst[i] = average[i] - (stdDev[i]/10);
		maxEst[i] = average[i] + (stdDev[i]/10);
	}
//time to output
//average
	FILE *out;
	char lname[strlen(name)+9];
	strcpy(lname,name);
	strcat(lname,"avg.graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	char nameEnd[12];
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"avg%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,average[i]);
	}
	fclose(out);
	printf("average g(r) saved to %s\n",lname);
//min
	strcpy(lname,name);
	strcat(lname,"min.graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"min%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,minEst[i]);
	}
	fclose(out);
	printf("estimated minimum g(r) saved to %s\n",lname);
//max
	strcpy(lname,name);
	strcat(lname,"max.graph");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoid overwriting old results
	j = 1;
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"max%d.graph",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	for(i=0;i<xMaxi*dr;i++){
		fprintf(out,"%f	%lf\n",(float)i/dr,maxEst[i]);
	}
	fclose(out);
	printf("estimated maximum g(r) saved to %s\n",lname);
	free(distances);
	return;
}
