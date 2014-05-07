#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>

double discretegeometricaverage(int T, int M, int S0, int K, double sigma, double r){
	double dt=(double)T/M;

	double T1=T-(M*(M-1)*(4*M+1))/(6.*M*M)*dt;
	double T2=T-(M-1)/2.*dt;

	double A=exp(-r*(T-T2)-sigma*sigma*(T2-T1)/2.);
	double d=(log((double)S0/K)+(r-0.5*sigma*sigma)*T2)/(sigma*sqrt(T1));

	return S0*A*gsl_cdf_gaussian_P(d+sigma*sqrt(T1),1)-K*exp(-r*T)*gsl_cdf_gaussian_P(d,1);
}

double continuousgeometricaverage(int T, int szero, int K, double sigma, double r){
	double d=(log((double)szero/K)+0.5*(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T/3.));

	return szero*exp(-0.5*(r+1./6*sigma*sigma)*T)*gsl_cdf_gaussian_P(d+sigma*sqrt(T/3.),1)-K*exp(-r*T)*gsl_cdf_gaussian_P(d,1);
}

double simulation(int T, int M, int S0, int K, double sigma, double mu, double r,gsl_rng* rng){

	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=S0;
	w[0]=gsl_ran_ugaussian(rng);
	double prod=1;

	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(i*dt-(i-1)*dt)*gsl_ran_ugaussian(rng);

		s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*dt+sigma*w[i]);
		prod*=s[i];
	}

	return std::max(pow(prod,1./M)-K,(double)0);
}

int main(){
	std::ofstream file;
	file.open("data.dat");

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	int n=100000;
	double sum=0;
	for(int i=0;i<n;i++)
		sum+=simulation(1,10,10,10,0.25,0,0.1,r);

	printf("%f\n%f\n",discretegeometricaverage(1,10,10,10,0.25,0.1),sum/n);

	file.close();


	return 1;
}