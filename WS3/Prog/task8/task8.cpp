#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <fstream>
#define _USE_MATH_DEFINES
#include <vector>
using std::vector;

double f(double x){
	return 1+0.1*exp(x/2.);
}

//tensor
double tensor(int level[], vector<vector<double> > nodes, vector<vector<double> > weights, int d, std::ofstream &file){
	//weights[i][j]: weight of j-th node on level i
	//nodes[i][j]: j-th node on level i 
	int k[d];

	for(int i=1;i<=d;i++)
		k[i]=1;

	double sum=0;
	double weight=1;
	double prod=1;

	while(1){
		//compute current weight and evaluate f
		for(int i=1;i<=d;i++){
			weight*=weights[level[i]][k[i]];
			prod*=f(nodes[level[i]][k[i]]);
		}
		sum+=weight*prod;

		//reset weight
		weight=1;
		prod=1;

		//output grid
		/*for(int i=1;i<=d;i++){
			file << nodes[level[i]][k[i]]<< " ";
		}
		file << "\n";*/

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

double sumall(int level, vector<vector<double> > nodes, vector<vector<double> > weights, int d){
	std::ofstream file;
	file.open("test");
	int k[d];

	double sum=0;

	for(int i=1;i<=d;i++)
		k[i]=1;

	while(1){
		sum+=tensor(k,nodes,weights,d,file);

		for(int j=1;j<=d;j++){
			k[j]++;

			//we want to get all possible combinations of nodes for the given level
			if(k[j]>level)
			{
				if(j==d)
					return sum;
				k[j]=1;
			}
			else
				break;
		}
	}
	file.close();
}

double sumsimplex(int level, vector<vector<double> > nodes, vector<vector<double> > weights, int d){
	std::ofstream file;
	file.open("test");

	//sum over simplex
	int k[d+1];
	int S=d;
	for(int i=1;i<=d;i++)
		k[i]=1;

	double sum=0;

	while(1){
		//call tensor function with given levels
		sum+=tensor(k,nodes,weights,d,file);

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
	file.close();
}

double SparseGridTrapezoidal(int level, int d){
	//init nodes,weights
	vector<vector<double> > nodes;
	nodes.resize(level+1);
	vector<vector<double> > weights;
	weights.resize(level+1);
	for(int i=1;i<=level;i++){
		nodes[i].resize(pow(2,level));
		weights[i].resize(pow(2,level));
	}

	//precalculate nodes,weights for levels <=level
	for(int l=1;l<=level;l++)
		for(int i=1;i<=pow(2,l)-1;i++){
			nodes[l][i]=i/(pow(2,l));
			if(i==1 || i==pow(2,l)-1)
				weights[l][i]=(3/2.)/(pow(2,l));
			else
				weights[l][i]=1./(pow(2,l));
		}

	for(int l=level;l>0;l--)
		for(int i=1;i<=pow(2,l)-1;i++){
			if(i%2==0 && l>1)
				weights[l][i]-=weights[l-1][i/2];
		}

	//sum over simplex
	printf("%f\n",sumsimplex(level,nodes,weights,d));
}

double ProductRuleTrapezoidal(int level, int d){
	//init nodes,weights
	vector<vector<double> > nodes;
	nodes.resize(level+1);
	vector<vector<double> > weights;
	weights.resize(level+1);
	for(int i=1;i<=level;i++){
		nodes[i].resize(pow(2,level));
		weights[i].resize(pow(2,level));
	}

	//precalculate nodes,weights for levels <=level
	for(int l=1;l<=level;l++)
		for(int i=1;i<=pow(2,l)-1;i++){
			nodes[l][i]=i/(pow(2,l));
			if(i==1 || i==pow(2,l)-1)
				weights[l][i]=(3/2.)/(pow(2,l));
			else
				weights[l][i]=1./(pow(2,l));

		}

	for(int l=level;l>0;l--)
		for(int i=1;i<=pow(2,l)-1;i++){
			if(i%2==0 && l>1)
				weights[l][i]-=weights[l-1][i/2];
		}

	//sum over simplex
	printf("%f\n",sumall(level,nodes,weights,d));
}

double SparseGridCC(int level, int d){
	//init nodes,weights
	vector<vector<double> > nodes;
	nodes.resize(level+1);
	vector<vector<double> > weights;
	weights.resize(level+1);
	for(int i=1;i<=level;i++){
		nodes[i].resize(pow(2,level));
		weights[i].resize(pow(2,level));
	}

	//precalculate nodes,weights for levels <=level
	for(int l=1;l<=level;l++)
		for(int i=1;i<=pow(2,l)-1;i++){
			nodes[l][i]=0.5*(1-std::cos(M_PI*i/pow(2,l)));
			
			double sum=0;
			for(int j=1;j<=pow(2,l)/2;j++)
			sum+=1./(2*j-1)*std::sin((2*j-1)*i*M_PI/pow(2,l));

			weights[l][i]=2./pow(2,l)*sin(i*M_PI/pow(2,l))*sum;
		}

	for(int l=level;l>0;l--)
		for(int i=1;i<=pow(2,l)-1;i++){
			if(i%2==0 && l>1)
				weights[l][i]-=weights[l-1][i/2];
		}

	//sum over simplex
	printf("%f\n",sumsimplex(level,nodes,weights,d));
}

double ProductRuleCC(int level, int d){
	//init nodes,weights
	vector<vector<double> > nodes;
	nodes.resize(level+1);
	vector<vector<double> > weights;
	weights.resize(level+1);
	for(int i=1;i<=level;i++){
		nodes[i].resize(pow(2,level));
		weights[i].resize(pow(2,level));
	}

	//precalculate nodes,weights for levels <=level
	for(int l=1;l<=level;l++)
		for(int i=1;i<=pow(2,l)-1;i++){
			nodes[l][i]=0.5*(1-std::cos(M_PI*i/pow(2,l)));
			
			double sum=0;
			for(int j=1;j<=pow(2,l)/2;j++)
			sum+=1./(2*j-1)*std::sin((2*j-1)*i*M_PI/pow(2,l));

			weights[l][i]=2./pow(2,l)*sin(i*M_PI/pow(2,l))*sum;
		}

	for(int l=level;l>0;l--)
		for(int i=1;i<=pow(2,l)-1;i++){
			if(i%2==0 && l>1)
				weights[l][i]-=weights[l-1][i/2];
		}

	//sum over simplex
	printf("%f\n",sumall(level,nodes,weights,d));
}

int main(){
	ProductRuleCC(5,3);

	return 1;
}