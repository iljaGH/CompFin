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
double tensor(int level[], vector<vector<double> > nodes, vector<vector<double> > weights, int d){
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
	int k[d];

	double sum=0;

	for(int i=1;i<=d;i++)
		k[i]=1;

	while(1){
		sum+=tensor(k,nodes,weights,d);

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
}

double sumsimplex(int level, vector<vector<double> > nodes, vector<vector<double> > weights, int d){
	//sum over simplex
	int k[d+1];
	int S=d;
	for(int i=1;i<=d;i++)
		k[i]=1;

	double sum=0;

	while(1){
		//call tensor function with given levels
		sum+=tensor(k,nodes,weights,d);

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
	return sumsimplex(level,nodes,weights,d);
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
	return sumall(level,nodes,weights,d);
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
	return sumsimplex(level,nodes,weights,d);
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
	return sumall(level,nodes,weights,d);
}

double MC(int level, int d){

	int n=pow(2,level)-1;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	double sum=0;
	double prod=1;

	for(long int i=0;i<n;i++){
		for(int j=0;j<d;j++)
			prod*=f(gsl_rng_uniform(r));
		sum+=prod;
		prod=1;
	}

	return sum/n;
}

int isprime(int n){
	for(int i=2;i<=sqrt(n);i++)
		if(n%i==0)
			return 0;

	return 1;
}

void primes(int primes[], int n){
	int k=0,p=0;
	while(k<n){
		p++;
		if(isprime(p)){
			primes[k]=p;
			k++;
		}
	}
}

void vandercorput(double vdc[],int n, int p, int epsilon=12){
	vdc[0]=0;
	for(int i=1;i<n;i++){
		double z=1-vdc[i-1];
		double v=1./p;
		while(z*pow(10,epsilon)<v*pow(10,epsilon)+1){
			v=v/p;
		}
		vdc[i]=vdc[i-1]+(p+1)*v-1;
	}
}

void halton(double** halton, int n, int d){
	int prime[d+3];
	primes(prime,d+3);
	double vdc[n+500];

	for(int i=0;i<d;i++)
	{
		vandercorput(vdc,n+500,prime[i+1]);
		for(int j=0;j<n;j++)
			halton[j][i]=vdc[j+500];
	}
}

double QMC(int level, int d){
	int n=pow(2,level)-1;

	double sum=0;
	double prod=1;

	double *halt[n];
	for(int i=0;i<n;i++)
		halt[i]=new double[d];

	halton(halt,n,d);

	for(int i=0;i<n;i++){
		for(int j=0;j<d;j++)
			prod*=f(halt[i][j]);

		sum+=prod;
		prod=1;
	}

	return sum/n;
}

int main(){
	std::ofstream file;
	int maxlevel=10;

	for(int d=1;d<=8;d*=2){
		switch(d){
			case 1: file.open("d1.dat");break;
			case 2: file.open("d2.dat");break;
			case 4: file.open("d4.dat");break;
			case 8: file.open("d8.dat");break;
		}

		double expected=pow(1.12974,d);
		file << "#level MC QMC SGCC PGCC SGT PGT d="<<d<<"\n";
		for(int l=1;l<10;l++){
			printf("d=%i level %i\n",d,l);
			//PRODUCT RULE DAUERT SO LANGE AAAAAAAAAAAAAAAAAH!!!!
			double prodcc=0,prodtrap=0;
			if(l<3){
				prodcc=std::abs(ProductRuleCC(l,d)-expected)/expected;
				prodtrap=std::abs(ProductRuleTrapezoidal(l,d)-expected)/expected;
			}
			file << l << " "<<std::abs(MC(l,d)-expected)/expected<< " "<<std::abs(QMC(l,d)-expected)/expected<< " "<<std::abs(SparseGridCC(l,d)-expected)/expected<< " "<<prodcc<< " "<<std::abs(SparseGridTrapezoidal(l,d)-expected)/expected<< " "<<prodtrap<< "\n";
		}

		file.close();
	}


	return 1;
}