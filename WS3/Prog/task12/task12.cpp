#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <fstream>
#define _USE_MATH_DEFINES
#include <vector>
#include <string>
using std::vector;

double f(double x){
	return 1+0.1*exp(x/2.);
}

//tensor
double tensor(int level[], int d){
	//weights[i][j]: weight of j-th node on level i
	//nodes[i][j]: j-th node on level i 
	int k[d];

	for(int i=1;i<=d;i++)
		k[i]=1;

	double sum=0;

	while(1){

		sum++;

		for(int j=1;j<=d;j++){
			k[j]++;

			//we want to get all possible combinations of nodes for the given levels
			if(k[j]>pow(2,level[j])-1)
			{
				if(j==d)
					return sum;
				k[j]=1;
			}
			else
				break;
		}
	}
}

double sumsimplex(int level, int d){
	//sum over simplex
	int k[d+1];
	int S=d;
	for(int i=1;i<=d;i++)
		k[i]=1;

	double sum=0;

	while(1){
		//call tensor function with given levels
		sum+=tensor(k,d);

		for(int j=1;j<=d;j++){
			k[j]++;
			S++;
			if(S>d+level-1){
				if(j==d)
					return sum;
				S=S-(k[j]-1);
				k[j]=1;
			} else{
				break;
			}

		}
	}
}

int main(){
	std::ofstream file;
	int level=4;
	file.open("data.dat");

	for(int d=1;d<=10;d++){
		printf("%i\n",d);
		file << d << " " << sumsimplex(level,d)<<" " << pow(pow(2,level)-1,d)<<"\n";

	}
	file.close();


	return 1;
}