#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>

int T=1;
int M=2;
int K=10;
int szero=10;
double r=0.1;
double sigma=0.25;

double randomwalk(double z[]){

	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=szero;
	w[0]=0;

	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(dt)*z[i-1];
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);
	}

	return std::max(0.5*(s[1]+s[2])-K,(double) 0);
}

double f(double x[]){
	for(int i=0;i<M;i++){
		x[i]=gsl_cdf_gaussian_Pinv(x[i],1);
	}


		return randomwalk(x);

}

double discretearithmeticaverage(int T, int M, int S0, int K, double sigma, double r,gsl_rng* rng){

	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=S0;
	w[0]=gsl_ran_ugaussian(rng);
	double sum=0;

	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(i*dt-(i-1)*dt)*gsl_ran_ugaussian(rng);

		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);
		sum+=s[i];
	}

	return std::max(sum/M-K,(double)0);
}

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
	file.open("task5.dat");
	int num=100;
	double delta=1./num;

	for(int i=1;i<num;i++)
		for(int j=1;j<num;j++){
			double x[M];
			x[0]=i*delta;
			x[1]=j*delta;

			file << delta*i << " " << delta*j << " "<<f(x)<<"\n";
		}


	file.close();


	return 1;
}
