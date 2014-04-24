#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <cmath>

void runTrial(gsl_rng* r,int N, int K, int szero, double mu, double sigma, int T, double deltat){

	int M=T/deltat;
	double w[M+1],s[M+1];

	double values[N];

	//get the call option value of every GBM.
	for(int j=0;j<N;j++){
		s[0]=szero;
		w[0]=gsl_ran_ugaussian(r);
		for(int i=1;i<=M;i++){
			w[i]=w[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
			
			s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w[i]);
		}
		values[j]=std::max(s[M]-K,(double)0);
	}

	//mean+variance
	double alpha=values[0];
	double beta=0;
	for(int i=1;i<N;i++){
		double gamma=values[i]-alpha;
		alpha+=gamma/(i+1);
		beta=beta+gamma*gamma*i/(i+1);
	}
	double esigma=sqrt(beta/(N-1));

	printf("deltat=%f,mean=%f,variance=%f\n",deltat,alpha,esigma);

}

int main(){
	int N=1000;
	int K=10;
	double mu=0.1;
	double sigma=0.2;
	int T=2;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	runTrial(r,N,K,10,mu,sigma,T,0.2);
	runTrial(r,N,K,10,mu,sigma,T,0.4);
	runTrial(r,N,K,10,mu,sigma,T,0.5);
	runTrial(r,N,K,10,mu,sigma,T,1.0);
	runTrial(r,N,K,10,mu,sigma,T,2.0);
}