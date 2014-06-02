#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <fstream>
#include <ctime>
#define _USE_MATH_DEFINES
#include <vector>
#include <string>
using std::vector;

//keine Lust in den Funktionen rumzuwühlen, deshalb globale Variablen
int T=1;
int M=64;
int K=10;
int szero=10;
double r=.1;
double sigma=.25;
int rndwlk=1;

double randomwalk(double z[]){

	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=szero;
	w[0]=0;
	double prod=1;

	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(dt)*z[i];
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);

		prod*=s[i];
	}

	return std::max(pow(prod,1./M)-K,(double)0);
}

double brownianbridge(double z[]){
	gsl_rng* rng;
	rng=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,time(NULL));

	double w[M+1],s[M+1];
	double dt=(double) T/M;

	s[0]=szero;
	w[0]=0;
	w[M]=sqrt(T)*z[0];
	double prod=1;

	//for each level do: if nodeindex is odd (not calculated yet): use formula
	int i=1;
	for(int l=1;pow(2,l)<=M;l++)
		for(int k=0;k*pow(2,-l)<T;k++)
			if(k%2==1)
				w[(int)(k*pow(2,-l)/dt)]=0.5*(w[(int)((k-1)*pow(2,-l)/dt)]+w[(int)((k+1)*pow(2,-l)/dt)])+sqrt(0.5*pow(2,-l)*T)*z[i++];

	for(int i=1;i<=M;i++){
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);
		prod*=s[i];
	}

	return std::max(pow(prod,1./M)-K,(double)0);
}

double f(double x[]){
	for(int i=0;i<M;i++)
		x[i]=gsl_cdf_gaussian_Pinv(x[i],1);

	if(rndwlk)
		return randomwalk(x);
	else
		return brownianbridge(x);
}

double tensor(int level[], vector<vector<double> > nodes, vector<vector<double> > weights, int d){
	//weights[i][j]: weight of j-th node on level i
	//nodes[i][j]: j-th node on level i 
	int k[d];

	for(int i=1;i<=d;i++)
		k[i]=1;

	double sum=0;
	double weight=1;

	while(1){
		//compute current weight and evaluate f
		double x[d];
		for(int i=1;i<=d;i++){
			weight*=weights[level[i]][k[i]];
			x[i-1]=nodes[level[i]][k[i]];
		}
		sum+=weight*f(x);

		//reset weight
		weight=1;

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
		double x[d];
		for(int j=0;j<d;j++)
			x[j]=gsl_rng_uniform(r);
		prod*=f(x);
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
		double x[d];
		for(int j=0;j<d;j++)
			x[j]=halt[i][j];
		prod*=f(x);

		sum+=prod;
		prod=1;
	}

	return sum/n;
}

double discretegeometricaverage(){
	double dt=(double)T/M;

	double T1=T-(M*(M-1)*(4*M+1))/(6.*M*M)*dt;
	double T2=T-(M-1)/2.*dt;

	double A=exp(-r*(T-T2)-sigma*sigma*(T2-T1)/2.);
	double d=(log((double)szero/K)+(r-0.5*sigma*sigma)*T2)/(sigma*sqrt(T1));

	return szero*A*gsl_cdf_gaussian_P(d+sigma*sqrt(T1),1)-K*exp(-r*T)*gsl_cdf_gaussian_P(d,1);
}


int main(){
	int maxlevel =10;
	std::ofstream file;
	file.open("mc.dat");
	double expected=discretegeometricaverage();
	file << "#level | error randomwalk | error brownian bridge\n";

	for(int i=1;i<maxlevel;i++){
		rndwlk=1;

		double res1=exp(-0.1)*MC(i,M);
		rndwlk=0;
		double res2=exp(-0.1)*MC(i,M);
		printf("%i %f %f %f\n",i,expected,res1,res2);
		file << i << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

	file.open("qmc.dat");
	file << "#level | error randomwalk | error brownian bridge\n";

	for(int i=1;i<maxlevel;i++){
		rndwlk=1;
		double res1=exp(-0.1)*QMC(i,M);
		rndwlk=0;
		double res2=exp(-0.1)*QMC(i,M);
		printf("%i %f %f %f\n",i,expected,res1,res2);
		file << i << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

	return 1;
}
