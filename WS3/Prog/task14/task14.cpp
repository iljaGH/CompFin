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

double randomwalk(int T, int M, int S0, int K, double sigma, double r){

	gsl_rng* rng;
	rng=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,time(NULL));

	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=S0;
	w[0]=gsl_ran_ugaussian(rng);

	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(dt)*gsl_ran_ugaussian(rng);
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);
	}
}

double brownianbridge(int T, int M, int S0, int K, double sigma, double r){
	gsl_rng* rng;
	rng=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,time(NULL));

	double w[M+1],s[M+1];
	double dt=(double) T/M;

	s[0]=S0;
	w[0]=gsl_ran_ugaussian(rng);
	w[M]=sqrt(T)*gsl_ran_ugaussian(rng);

	//for each level do: if nodeindex is odd (not calculated yet): use formula
	for(int l=1;pow(2,l)<=M;l++)
		for(int k=0;k*pow(2,-l)<T;k++)
			if(k%2==1)
				w[(int)(k*pow(2,-l)/dt)]=0.5*(w[(int)((k-1)*pow(2,-l)/dt)]+w[(int)((k+1)*pow(2,-l)/dt)])+sqrt(0.5*pow(2,-l)*T)*gsl_ran_ugaussian(rng);

	for(int i=1;i<=M;i++)
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);
}

int main(){
	brownianbridge(1,16,10,2,3,4);

	return 1;
}