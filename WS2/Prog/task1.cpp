#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <cmath>

int main(){
	int N=1000;
	int K=10;
	double mu=0.1;
	int T=2;
	double deltat=0.2;

	int M=T/deltat;
	double w[M+1],s[M+1];
	std::ofstream file;

	file.open("task1.dat");

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	double sum=0;
	for(double sigma=0;sigma<=0.8;sigma+=0.2){
		for(int j=1;j<=N;j++){
			s[0]=10;
			w[0]=gsl_ran_ugaussian(r);
			for(int i=1;i<=M;i++){
				w[i]=w[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
				
				s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w[i]);
			}
			sum+=std::max(s[M]-K,(double)0);
		}
		file << sigma <<" " <<sum/N<<"\n";
	}
	file.close();
}