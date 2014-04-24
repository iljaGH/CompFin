#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <cmath>

void runTrial(int T, double mu, double sigma, double deltat, gsl_rng* r,int n){
	std::ofstream file;

	if(n==1)
		file.open("task10_t1.dat");
	if(n==2)
		file.open("task10_t2.dat");

	int M=T/deltat;
	double w1[M+1], w2[M+1], w3[M+1],s1[M+1], s2[M+1], s3[M+1];

	s1[0]=10;
	s2[0]=10;
	s3[0]=10;
	w1[0]=gsl_ran_ugaussian(r);
	w2[0]=gsl_ran_ugaussian(r);
	w3[0]=gsl_ran_ugaussian(r);


	file << "#deltat w1 w2 w3 s1 s2 s3\n";
	file << "0 " << w1[0]<< " " << w2[0] << " " << w3[0] << " " << s1[0] <<" " << s2[0] <<" " <<s3[0] <<" "<< "\n";

	for(int i=1;i<=M;i++){
		w1[i]=w1[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
		w2[i]=w2[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
		w3[i]=w3[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
		s1[i]=s1[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w1[i]);
		s2[i]=s2[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w2[i]);
		s3[i]=s3[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w3[i]);
		file << i*deltat << " " <<w1[i]<<" " <<w2[i]<<" " <<w3[i]<< " " << s1[i] <<" " << s2[i]<<" " << s3[i]<< "\n";
	}

	file.close();

}

int main(){
	int T=2;
	double mu=0.1, sigma=0.2,deltat1=0.5, deltat2=0.01;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);

	runTrial(T,mu,sigma,deltat1,r,1);
	runTrial(T,mu,sigma,deltat2,r,2);


	return 0;
}
