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



int main(){
	std::ofstream file;
	file.open("data.dat");

	int T=1,S0=10,K=10;
	double sigma=.25,r=.1;

	double expected=continuousgeometricaverage(T,S0,K,sigma,r);
	file << "#M error\n";

	for(int i=1;i<16;i++)
		file << pow(2,i)<<" "<<std::abs(expected-discretegeometricaverage(T,pow(2,i),S0,K,sigma,r))/expected<<"\n";

	file.close();


	return 1;
}