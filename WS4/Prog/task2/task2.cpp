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
int M=128;
int K=10;
int szero=10;
double r=.05;
double sigma=.2;
int rndwlk=1;
int B=9;

double randomwalk(double z[]){

	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=szero;
	w[0]=0;

	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(dt)*z[i-1];
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);

		//down out
		if(s[i]<=B)
			return 0;
	}

	return std::max(s[M]-K,(double) 0);
}

double brownianbridge(double z[]){
	gsl_rng* rng;
	rng=gsl_rng_alloc(gsl_rng_mt19937);
	//gsl_rng_set(rng,time(NULL));

	double w[M+1],s[M+1];
	double dt=(double) T/M;

	s[0]=szero;
	w[0]=0;
	w[M]=sqrt(T)*z[0];

	//for each level do: if nodeindex is odd (not calculated yet): use formula
	int i=1;
	for(int l=1;pow(2,l)<=M;l++)
		for(int k=0;k*pow(2,-l)<T;k++)
			if(k%2==1)
				w[(int)(k*pow(2,-l)/dt)]=0.5*(w[(int)((k-1)*pow(2,-l)/dt)]+w[(int)((k+1)*pow(2,-l)/dt)])+sqrt(0.5*pow(2,-l)*T)*z[i++];

	for(int i=1;i<=M;i++){
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);
		//Down out
		if(s[i]<=B)
			return 0;
	}

	return std::max(s[M]-K,(double) 0);
}

double f(double x[]){
	for(int i=0;i<M;i++)
		x[i]=gsl_cdf_gaussian_Pinv(x[i],1);

	if(rndwlk)
		return randomwalk(x);
	else
		return brownianbridge(x);
}

double MC(int level, int d){

	int n=pow(2,level)-1;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	//gsl_rng_set(r,time(NULL));

	double sum=0;

	for(long int i=0;i<n;i++){
		double x[d];
		for(int j=0;j<d;j++)
			x[j]=gsl_rng_uniform(r);
		sum+=f(x);
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

	double *halt[n];
	for(int i=0;i<n;i++)
		halt[i]=new double[d];

	halton(halt,n,d);
	for(int i=0;i<n;i++){
		double x[d];
		for(int j=0;j<d;j++)
			x[j]=halt[i][j];
		double res=f(x);
		if(res!=res){}
			else
		sum=(sum*i+res)/(i+1);
	}

	return sum;
}

int main(){
	int maxlevel =17;
	std::ofstream file;

	file.open("mc_high.dat");
	file << "#nodes | error randomwalk | error brownian bridge\n";

	//double expected=exp(-r)*MC(20,M);
	//printf("%f\n",expected);
	double expected=0.89799;


	for(int i=1;i<=maxlevel;i++){
		rndwlk=1;
		double res1=exp(-r)*MC(i,M);
		rndwlk=0;
		double res2=exp(-r)*MC(i,M);
		printf("%i %f %f\n",i,res1,res2);
		file << pow(2,i)-1 << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

	file.open("qmc.dat");
	file << "nodes | error randomwalk | error brownian bridge\n";

	for(int i=1;i<=maxlevel;i++){
		rndwlk=1;
		double res1=exp(-r)*QMC(i,M);
		rndwlk=0;
		double res2=exp(-r)*QMC(i,M);
		printf("%i %f %f\n",i,res1,res2);
		if(res1!=res1 || res2!=res2)
			break;
		file << pow(2,i)-1 << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

	file.open("mc_low.dat");
	file << "#nodes | error randomwalk | error brownian bridge\n";

	expected=exp(-r)*MC(6,M);

	for(int i=1;i<=maxlevel;i++){
		rndwlk=1;
		double res1=exp(-r)*MC(i,M);
		rndwlk=0;
		double res2=exp(-r)*MC(i,M);
		printf("%i %f %f\n",i,res1,res2);
		file << pow(2,i)-1 << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

	return 1;
}
