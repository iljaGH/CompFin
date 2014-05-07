#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <fstream>
#define _USE_MATH_DEFINES

double f(double x){
	return 1;
}

double tensor(int level[], double** nodes, int d, std::ofstream &file){
	//weights[i][j]: weight of j-th node on level i
	//nodes[i][j]: j-th node on level i
	int k[d];

	for(int i=0;i<d;i++)
		k[i]=1;

	while(1){

		double prod=1;
		for(int i=0;i<d;i++){
			printf("%i ",k[i]);
			file << nodes[level[i]-1][k[i]-1] << " ";
		}
		for(int i=0;i<d;i++){
			printf(" %f ",nodes[level[i]-1][k[i]-1]);
			//file << nodes[level[i]][k[i]] << " ";
		}
		printf("\n");
		file << "\n";
		
		for(int j=0;j<d;j++){
			k[j]++;
			if(k[j]>pow(2,level[j])-1)
			{
				if(j==d-1){
					printf("------\n");
					return 1;
				}
				k[j]=1;
			}
			else
				break;
		}
	}
}

double SparseGridTrapezoidal(int level, int d){

	std::ofstream file;
	file.open("test");
	double *nodes[level];
	for(int i=0;i<level;i++)
		nodes[i]=new double[(int)pow(2,level)-1];

	for(int i=0;i<level;i++)
		for(int n=0;n<pow(2,i+1)-1;n++)
		{
			if(i==2)
				printf("%f ",(n+1)/(double)(pow(2,i+1)));
			nodes[i][n]=(n+1)/(double)(pow(2,i+1));
		}
		printf("\n");

	int k[d];
	int S=d;
	for(int i=0;i<d;i++)
		k[i]=1;

	while(1){

		tensor(k,nodes,d,file);

		for(int j=0;j<d;j++){
			k[j]++;
			S++;
			if(S>d+level-1){
				if(j==d-1)
					return 1;
				S=S-(k[j]-1);
				k[j]=1;
			} else{
				break;
			}

		}
	}
	file.close();
}

int main(){
	SparseGridTrapezoidal(5,2);

	return 1;
}