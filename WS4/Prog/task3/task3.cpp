#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>
#include <gsl/gsl_cdf.h>
#define _USE_MATH_DEFINES

double szero=10;
int K=10;
int T=1;
double sigma=0.2;
double r=0.05;
double B=9;

double d(double p, double q){
	return (log(p/q)+(r-sigma*sigma*0.5)*T)/(sigma*sqrt(T));
}

double v_bs(double p, double q){
	return p*gsl_cdf_ugaussian_P(d(p,q)+sigma*sqrt(T))-K*exp(-r*T)*gsl_cdf_ugaussian_P(d(p,q));
}

double fairprice_downout(){
	double Bhat=std::max(B,(double)K);
	double Z = pow((B/szero),(2*r)/(sigma*sigma)-1);

	return v_bs(szero,Bhat)-Z*v_bs(B*B/szero,Bhat)+(Bhat-K)*exp(-r*T)*(gsl_cdf_ugaussian_P(d(szero,Bhat))-Z*gsl_cdf_ugaussian_P(d(B*B/szero,Bhat)));
}

int main(){
	std::ofstream file;
	file.open("task3.dat");
	for(B=0;B<10;B+=0.1){
		printf("%f %f\n",B,fairprice_downout());
		file << B << " " << fairprice_downout() << "\n";
	}
	file.close();

	return 1;
}
