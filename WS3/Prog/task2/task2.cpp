#include <gsl/gsl_cdf.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

double discretegeometricaverage(int T, int M, int szero, int K, double sigma, double r){
	double dt=(double)T/M;

	double T1=T-(M*(M-1)*(4*M+1))/(6*M*M)*dt;
	double T2=T-(M-1)/2*dt;

	double A=exp(-r*(T-T2)-sigma*sigma*(T2-T1)/2);
	double d=(log(szero/K)+(r-0.5*sigma*sigma)*T2)/(sigma*sqrt(T1));

	return szero*A*gsl_cdf_gaussian_P(d+sigma*sqrt(T1),1)-K*exp(-r*T)*gsl_cdf_gaussian_P(d,1);
}

double continuousgeometricaverage(int T, int szero, int K, double sigma, double r){
	double d=(log(szero/K)+0.5*(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T/3));

	return szero*exp(-0.5*(r+(double)1/6*sigma*sigma)*T)*gsl_cdf_gaussian_P(d+sigma*sqrt(T/3),1)-K*exp(-r*T)*gsl_cdf_gaussian_P(d,1);
}

int main(){
	printf("%f",gsl_cdf_gaussian_P(0.3,1));

	return 1;
}