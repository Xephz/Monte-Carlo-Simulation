#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include"coreMCfunc.h"
#include"particleIO.h"

void Populate_Random(particle *particles,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles){
	int i;
	for(i=0;i<numParticles;i++){
		(particles+i)->x = PRand(randoms,currentRandom)*xMax;
		(particles+i)->y = PRand(randoms,currentRandom)*yMax;
		(particles+i)->z = PRand(randoms,currentRandom)*zMax;
		(particles+i)->cluster = 0;
	}
	return;
}

void Populate_Cubic(particle *particles,double xMax,int numParticles){
	int i = 0,x,y,z;
	int numIntervals = cbrt(numParticles);
	double intervalWidth;
	intervalWidth = (xMax/numIntervals);
	for(x=0;x<numIntervals;x++){
		for(y=0;y<numIntervals;y++){
			for(z=0;z<numIntervals;z++){
				(particles+i)->x = x*intervalWidth;
				(particles+i)->y = y*intervalWidth;
				(particles+i)->z = z*intervalWidth;	
				(particles+i)->cluster = 0;
				i++;
			}
		}
	}
	return;
}

void Populate_Cluster_fcc(particle *particles, double a, double spread,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles){
	int i = 0, xCount, yCount, zCount;
	int numXMoves = floor(xMax/a);
	int numYMoves = floor(yMax/a);
	int numZMoves = floor(zMax/a);
	double displaceDist;//distance to displace each particle from lattice site
	double displaceAngle1;//angle the displacement makes in x-y plane with x axis
	double displaceAngle2;//angle the displacement makes in x-z plane with z axis
	while (i<numParticles){
		//these 3 loops navigate to each lattice point
		for(xCount=0;xCount<numXMoves;xCount++){
			for(yCount=0;yCount<numYMoves;yCount++){
				for(zCount=0;zCount<numZMoves;zCount++){
					//add the corner particle
					displaceDist = (PRand(randoms,currentRandom)-0.5)*2*spread;
					displaceAngle1 = PRand(randoms,currentRandom)*M_PI;
					displaceAngle2 = PRand(randoms,currentRandom)*M_PI;
					(particles+i)->x = a*xCount+displaceDist*cos(displaceAngle1)*sin(displaceAngle2);
					(particles+i)->y = a*yCount+displaceDist*sin(displaceAngle1);
					(particles+i)->z = a*zCount+displaceDist*cos(displaceAngle2);
					i++;
					//add the x=0 particle
					displaceDist = (PRand(randoms,currentRandom)-0.5)*2*spread;
					displaceAngle1 = PRand(randoms,currentRandom)*M_PI;
					displaceAngle2 = PRand(randoms,currentRandom)*M_PI;
					(particles+i)->x = a*xCount+displaceDist*cos(displaceAngle1)*sin(displaceAngle2);
					(particles+i)->y = a*(yCount+0.5)+displaceDist*sin(displaceAngle1);
					(particles+i)->z = a*(zCount+0.5)+displaceDist*cos(displaceAngle2);
					i++;
					//add the y=0 particle
					displaceDist = (PRand(randoms,currentRandom)-0.5)*2*spread;
					displaceAngle1 = PRand(randoms,currentRandom)*M_PI;
					displaceAngle2 = PRand(randoms,currentRandom)*M_PI;
					(particles+i)->x = a*(xCount+0.5)+displaceDist*cos(displaceAngle1)*sin(displaceAngle2);
					(particles+i)->y = a*yCount+displaceDist*sin(displaceAngle1);
					(particles+i)->z = a*(zCount+0.5)+displaceDist*cos(displaceAngle2);
					i++;
					//add the z=0 particle
					displaceDist = (PRand(randoms,currentRandom)-0.5)*2*spread;
					displaceAngle1 = PRand(randoms,currentRandom)*M_PI;
					displaceAngle2 = PRand(randoms,currentRandom)*M_PI;
					(particles+i)->x = a*(xCount+0.5)+displaceDist*cos(displaceAngle1)*sin(displaceAngle2);
					(particles+i)->y = a*(yCount+0.5)+displaceDist*sin(displaceAngle1);
					(particles+i)->z = a*zCount+displaceDist*cos(displaceAngle2);
					i++;
				}
			}
		}
	}
	return;
}

void Populate_Cluster_bcc(particle *particles, double a, double spread,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles){
	int i = 0, xCount, yCount, zCount;
	int numXMoves = floor(xMax/a);
	int numYMoves = floor(yMax/a);
	int numZMoves = floor(zMax/a);
	double displaceDist;//distance to displace each particle from lattice site
	double displaceAngle1;//angle the displacement makes in x-y plane with x axis
	double displaceAngle2;//angle the displacement makes in x-z plane with z axis
	while (i<numParticles){
		//these 3 loops navigate to each lattice point
		for(xCount=0;xCount<numXMoves;xCount++){
			for(yCount=0;yCount<numYMoves;yCount++){
				for(zCount=0;zCount<numZMoves;zCount++){
					//add the corner particle 
					displaceDist = (PRand(randoms,currentRandom)-0.5)*2*spread;
					displaceAngle1 = PRand(randoms,currentRandom)*M_PI;
					displaceAngle2 = PRand(randoms,currentRandom)*M_PI;
					(particles+i)->x = a*xCount+displaceDist*cos(displaceAngle1)*sin(displaceAngle2);
					(particles+i)->y = a*yCount+displaceDist*sin(displaceAngle1);
					(particles+i)->z = a*zCount+displaceDist*cos(displaceAngle2);
					i++;
					//add the centre particle
					displaceDist = (PRand(randoms,currentRandom)-0.5)*2*spread;
					displaceAngle1 = PRand(randoms,currentRandom)*M_PI;
					displaceAngle2 = PRand(randoms,currentRandom)*M_PI;
					(particles+i)->x = a*(xCount+0.5)+displaceDist*cos(displaceAngle1)*sin(displaceAngle2);
					(particles+i)->y = a*(yCount+0.5)+displaceDist*sin(displaceAngle1);
					(particles+i)->z = a*(zCount+0.5)+displaceDist*cos(displaceAngle2);
					i++;
				}
			}
		}
	}
	return;
}

void Output_print(particle *particles,int numToOutput){
	int i;
	for(i=0;i<numToOutput;i++){
		printf("x=%lf y=%lf z=%lf\n",(particles+i)->x,(particles+i)->y,(particles+i)->z);
	}
	return;
}

void Output_jmol(particle *particles,char *name,double scaling,int numToOutput,double xMax,double yMax,double zMax){
	FILE *out;
	char lname[strlen(name)+4];
	strcpy(lname,name);
	strcat(lname,".xyz");
	//this section checks if the file already exists and if so renames untill an unused filename is found before opening, this avoids overwriting old results
	int j = 1;
	char nameEnd[7];
	while( access( lname, F_OK ) != -1){
		sprintf(nameEnd,"%d.xyz",j);
		strcpy(lname,name);
		strcat(lname,nameEnd);
		j++;
		//as we have only allowed 3 digits for j in nameEnd, this will prevent a seg fault if we are writing copious amounts of data
		if(j>999){
			break;
		}
	}
	out=fopen(lname,"w");
	fprintf(out,"%d\n\n",numToOutput+4);//adds header
	fprintf(out,"H 0.000000 0.000000 0.000000\nH 0.000000 %lf 0.000000\nH %lf 0.000000 0.000000\nH 0.000000 0.000000 %lf\n",yMax*scaling,xMax*scaling,zMax*scaling);//adds axis markers
	int i;
	for(i=0;i<numToOutput;i++){
		switch(((particles+i)->cluster)%14){//different elements are randomly assigned to each cluster as they will then appear as different colours in jmol and be easier to distinguish
			case 0:
			fprintf(out,"Fe %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 1:
			fprintf(out,"O %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 2:
			fprintf(out,"N %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 3:
			fprintf(out,"C %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 4:
			fprintf(out,"S %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 5:
			fprintf(out,"Cl %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 6:
			fprintf(out,"Cu %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 7:
			fprintf(out,"Kr %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 8:
			fprintf(out,"Br %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 9:
			fprintf(out,"Po %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 10:
			fprintf(out,"Rf %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 11:
			fprintf(out,"Cs %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 12:
			fprintf(out,"B %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 13:
			fprintf(out,"F %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
			case 14:
			fprintf(out,"Ra %lf %lf %lf\n",(particles+i)->x*scaling,(particles+i)->y*scaling,(particles+i)->z*scaling);
			continue;
		}
		printf("write error:output is probably wrong\n");
	}
	fclose(out);
	printf("%d particles saved to %s\n",numToOutput,lname);
	return;
}

void Input_jmol(particle *particles,char *name,double scaling,double *randoms,long int *currentRandom,double xMax,double yMax,double zMax,int numParticles){
	FILE *in;
	int error = 0;
	char lname[strlen(name)+4];
	strcpy(lname,name);
	strcat(lname,".xyz");
	if(access( lname, F_OK ) == -1){
		printf("file not found, populating randomly\n");
		Populate_Random(particles,randoms,currentRandom,xMax,yMax,zMax,numParticles);
		return;
	}
	in=fopen(lname,"r");
	int file_numParticles;
	if(fscanf(in,"%d\n\n",&file_numParticles)==EOF){
		printf("read error: can't read number of particles, check formatting\n");
	}
	file_numParticles -= 4;
	if(file_numParticles!=numParticles){
		printf("read error: number of particles inconsistant\n");
		error = 1;
		return;
	}
	//reads & sets axis maximums
	if(fscanf(in,"H 0.000000 0.000000 0.000000\nH 0.000000 %lf 0.000000\n",&yMax)==EOF){
		printf("read error: can't read yMax, check formatting\n");
		error = 1;
	}
	else{//correct for scaling on successful read
		yMax = yMax/scaling;
	}
	if(fscanf(in,"H %lf 0.000000 0.000000\n",&xMax)==EOF){
		printf("read error: can't read xMax, check formatting\n");
		error = 1;
	}
	else{//correct for scaling on successful read
		xMax = xMax/scaling;
	}
	if(fscanf(in,"H 0.000000 0.000000 %lf\n",&zMax)==EOF){
		printf("read error: can't read zMax, check formatting\n");
		error = 1;
	}
	else{//correct for scaling on successful read
		zMax = zMax/scaling;
	}
	//reads values
	int i;
	for(i=0;i<numParticles;i++){
		//this loop is to get past the type of particle which will be 1 or 2 characters long followed by a space
		char check = '\0';
		while(check != ' '){
			if(fscanf(in,"%c",&check)==EOF){
				error = 1;
			}
		}
		if(fscanf(in," %lf %lf %lf\n",&(particles+i)->x,&(particles+i)->y,&(particles+i)->z)==EOF){
			printf("read error: can't read particles, check formatting\n");
			error = 1;
		}
		else{//correct for scaling on successful read
			(particles+i)->x = (particles+i)->x/scaling;
			(particles+i)->y = (particles+i)->y/scaling;
			(particles+i)->z = (particles+i)->z/scaling;
		}
	}
	if(error){
		printf("Read was likely unsuccessful\n");
	}
	else{
		printf("%d particles read from %s successfully\n",numParticles,lname);
	}
	return;
}
