#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <fstream>
#define _USE_MATH_DEFINES
#include <vector>
using std::vector;

double f(double x){
	return 1;
}

//tensor
double tensor(int level[], int d,vector<vector<int> > &list){
	//weights[i][j]: weight of j-th node on level i shifted down by 1
	//nodes[i][j]: j-th node on level i shifted down by 1
	int k[d];

	for(int i=0;i<d;i++)
		k[i]=1;

	while(1){
		for(int j=0;j<d;j++)
			list[]

		//j is dimension shifted down by 1
		for(int j=0;j<d;j++){
			k[j]++;

			//we want to get all possible combinations of nodes for the given levels
			if(k[j]>pow(2,level[j])-1)
			{
				if(j==d-1)
					return 1;
				k[j]=1;
			}
			else
				break;
		}
	}
}

double sumsimplex(int level, int d){
	std::ofstream file;
	file.open("test");
	vector<vector<int> > list;

	list.resize(pow(2,level)-1);
	for(int i=0;i<pow(2,level)-1;i++)
		list[i].resize(pow(2,level)-1);

	//sum over simplex
	int k[d];
	int S=d;
	for(int i=0;i<d;i++)
		k[i]=1;

	while(1){
		//call tensor product formula
		tensor(k,d,list);

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
	sumsimplex(7,3);

	return 1;
}