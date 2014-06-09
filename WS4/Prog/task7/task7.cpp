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
double T=0.08767;
int K=50;
double szero=51.25;
double r=.05;
double sigma=.2;
double V=2;


double dphi(double x){
	return 1./sqrt(2*M_PI)*exp(-x*x*0.5);
}

double greek(double sigma1){
	return szero*sqrt(T)*dphi((log((double)szero/K)+(r+0.5*sigma1*sigma1)*T)/(sigma1*sqrt(T)));
}


double v(double sigma1){
	double d1=(log(szero/K)+(r+sigma1*sigma1*0.5)*T)/(sigma1*sqrt(T));
	double d2=d1-sigma1*sqrt(T);
	return szero*gsl_cdf_ugaussian_P(d1)-K*exp(-r*T)*gsl_cdf_ugaussian_P(d2);
}

double impliedsigma(){
	double impl=2*V/(sqrt(T)*szero);

	for(int i=0;i<300;i++)
		impl=impl-(v(impl)-V)/greek(impl);
	return impl;
}

int main(){
	int maxlevel =16;

	printf("%f\n",impliedsigma());

	printf("%f\n",v(impliedsigma()));

	return 1;
}
