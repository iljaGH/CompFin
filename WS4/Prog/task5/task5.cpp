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
double r=0.05;
double sigma=0.2;

double randomwalk(double z[]){
	double w[M+1],s[M+1];
	double dt=(double)T/M;

	s[0]=szero;
	w[0]=0;

	double max=0;
	for(int i=1;i<=M;i++){
		w[i]=w[i-1]+sqrt(dt)*gsl_cdf_ugaussian_Pinv(z[i-1]);
		s[i]=s[0]*exp((r-0.5*sigma*sigma)*i*dt+sigma*w[i]);

		if(s[i]>max)
			max=s[i];
	}

	return std::max(max-K,(double) 0);
}

double f(double x[]){
	return randomwalk(x);

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
