#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <cmath>
#include <time.h>

void runTrial(gsl_rng* r,int N, int K, int szero, double mu, double sigma, int T, double deltat, int trial){

	std::ofstream file;

	switch(trial){
		case 1:file.open("task4_1.dat");break;
		case 2:file.open("task4_2.dat");break;
		case 3:file.open("task4_3.dat");break;
		case 4:file.open("task4_4.dat");break;
		case 5:file.open("task4_5.dat");break;
	}

	int M=T/deltat;
	double w[M+1],s[M+1];

	double values[N];

	//get the call option value of every GBM
	for(int j=0;j<N;j++){
		s[0]=szero;
		w[0]=gsl_ran_ugaussian(r);
		for(int i=1;i<=M;i++){
			w[i]=w[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
			
			s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w[i]);
		}
		values[j]=std::max(s[M]-K,(double)0);
	}
	//mean
	double alpha=values[0];
	file << "1 " << alpha <<"\n";
	for(int i=1;i<N;i++){
		double gamma=values[i]-alpha;
		alpha+=gamma/(i+1);
		file << i+1 <<" " << alpha <<"\n";
	}

	file.close();
}


int main(){
	double szero=10;
	int N=1000000;
	int K=10;
	double mu=0.1;
	double sigma=0.2;
	int T=1;
	double deltat=1;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	for(int i=1;i<=5;i++)
		runTrial(r,N,K,szero,mu,sigma,T,deltat,i);
}