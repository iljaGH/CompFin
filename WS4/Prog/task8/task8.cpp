#include <cstdio>
#include <iostream>
#include <iomanip>
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
double T=1./365;
double K=10050;
double szero=9989;
double r=.03;
double V=0.07;

double dphi(double x){
	return 1./sqrt(2*M_PI)*exp(-(x*x)*0.5);
}

double vega(double sigma){
	double d1=(log(szero/K)+(r+sigma*sigma*0.5)*T)/(sigma*sqrt(T));

	return szero*sqrt(T)*dphi(d1);
}


double BlackScholes(double sigma){
	double d1=(log(szero/K)+(r+sigma*sigma*0.5)*T)/(sigma*sqrt(T));
	double d2=d1-sigma*sqrt(T);
	return szero*gsl_cdf_ugaussian_P(d1)-K*exp(-r*T)*gsl_cdf_ugaussian_P(d2);
}

double ImpliedVolNewton(double start=(V/szero) / (0.398 * sqrt(T))) {
	double impl=start;
	for(int i=0;i<20;i++)
		impl=impl-(BlackScholes(impl)-V)/vega(impl);

	return impl;
}

double ImpliedVolBisection(double tol=0.0000000001){
	double a=0,b=2;
	int n=0;
	while(1){
		double mid=a+0.5*(b-a),	Bmid=BlackScholes(a+0.5*(b-a))-V;
		if(Bmid<0)
			a=mid;
		else
			b=mid;
		n++;
		if(n==20)
			break;

		if(std::abs(Bmid)<=tol)
			return mid;
	}
}

int main(){
	std::ifstream infile;
	std::ofstream outfile;
	string dump;
	infile.open("data.dat");
	getline(infile,dump);
	infile >> T >> szero;

	outfile.open("out.dat");
	while(infile >> K >> V){
		cout << K << " "<<ImpliedVolBisection()<<"\n";
		outfile << K << " " <<ImpliedVolBisection()<<"\n";
	}

	outfile.close();
	infile.close();

	return 1;
}
