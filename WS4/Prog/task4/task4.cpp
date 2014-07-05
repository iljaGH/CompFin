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
using namespace std;

//keine Lust in den Funktionen rumzuwühlen, deshalb globale Variablen
int T=1;
int M=128;
int K=10;
int szero=10;
double r=0.05;
double sigma=0.2;
int rndwlk=1;
double B=9;

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

	gsl_rng* rng;
	rng=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,time(NULL));

	double sum=0;

	for(long int i=0;i<n;i++){
		double x[d];
		for(int j=0;j<d;j++)
			x[j]=gsl_rng_uniform(rng);
		sum+=f(x);
	}


	return sum/n;
}

double d(double p, double q){
	return (log(p/q)+(r-sigma*sigma*0.5)*T)/(sigma*sqrt(T));
}

double v_bs(double p, double q){
	return p*gsl_cdf_ugaussian_P(d(p,q)+sigma*sqrt(T))-K*exp(-r*T)*gsl_cdf_ugaussian_P(d(p,q));
}

double BlackScholes(double szero, double K){
	double d1=(log(szero/K)+(r+sigma*sigma*0.5)*T)/(sigma*sqrt(T));
	double d2=d1-sigma*sqrt(T);
	return szero*gsl_cdf_ugaussian_P(d1)-K*exp(-r*T)*gsl_cdf_ugaussian_P(d2);
}

double fairprice_downout(){
	double Bhat=std::max((double)B,(double)K);
	double Z = pow(((double)B/szero),(2*r)/(sigma*sigma)-1);

	return BlackScholes(szero,Bhat)-Z*BlackScholes(B*B/szero,Bhat)+(Bhat-K)*exp(-r*T)*(gsl_cdf_ugaussian_P(d(szero,Bhat))-Z*gsl_cdf_ugaussian_P(d(B*B/szero,Bhat)));
}

int main(){
	int maxlevel =16;
	std::ofstream file;

	double expected=fairprice_downout();
	printf("%f\n",expected);
	M=4;

	file.open("mc_4.dat");
	file << "#nodes | error randomwalk | error brownian bridge\n";

	for(int i=1;i<=maxlevel;i++){
		rndwlk=1;
		double res1=exp(-r)*MC(i,M);
		rndwlk=0;
		double res2=exp(-r)*MC(i,M);
		printf("%i %f %f\n",i,res1,res2);
		file << pow(2,i)-1 << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

		M=64;

	file.open("mc_64.dat");
	file << "#nodes | error randomwalk | error brownian bridge\n";

	for(int i=1;i<=maxlevel;i++){
		rndwlk=1;
		double res1=exp(-r)*MC(i,M);
		rndwlk=0;
		double res2=exp(-r)*MC(i,M);
		printf("%i %f %f\n",i,res1,res2);
		file << pow(2,i)-1 << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();


		M=256;

	file.open("mc_256.dat");
	file << "#nodes | error randomwalk | error brownian bridge\n";

	for(int i=1;i<=maxlevel;i++){
		rndwlk=1;
		double res1=exp(-r)*MC(i,M);
		rndwlk=0;
		double res2=exp(-r)*MC(i,M);
		printf("%i %f %f\n",i,res1,res2);
		file << pow(2,i)-1 << " "<< std::abs(expected-res1)/expected << " "<< std::abs(expected-res2)/expected << "\n";
	}

	file.close();

		M=1024;

	file.open("mc_1024.dat");
	file << "#nodes | error randomwalk | error brownian bridge\n";

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
