#include <cstdio>
#include <cstdlib>
#include <cmath>

//integrand?!
double f(double x){
	return 0;
}

double tensor(int level, int n, double** weights, double** nodes, int d){
	//weights[i][j]: weight of j-th node on level i
	//nodes[i][j]: j-th node on level i

	int k[d];

	for(int i=0;i<d;i++)
		k[i]=1;

	double sum=0;

	while(1){

		double weight=1,prod=1;
		for(int i=0;i<d;i++){
			weight*=weight[level][k[i]];
			prod*=f(nodes[level][k[i]]);
		}
		sum+=weight*prod;
		
		for(int j=0;j<d;j++){
			k[j]++;
			if(k[j]>n)
			{
				if(j==d-1)
					return sum;
				k[j]=1;
			}
			else
				break;
		}
	}
}

int main(){
	int level=3;
	int n=pow(2,level)-1;
	double weights[n];
	double nodes[n];
	int d=3;

	tensor(level,n,d);

	return 1;
}
