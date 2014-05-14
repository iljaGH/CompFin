#include "includes.hpp"

void runTrial(gsl_rng* r,int N, int K, int szero, double mu, double sigma, int T, double deltat, int trial){

	std::ofstream file;

	switch(trial){
		case 1:file.open("task4_1.dat");break;
		case 2:file.open("task4_2.dat");break;
		case 3:file.open("task4_3.dat");break;
		case 4:file.open("task4_4.dat");break;
		case 5:file.open("task4_5.dat");break;
	}

	int M=T/deltat;
	double w[M+1],s[M+1];

	double values = 0.0;
	double alpha = values;

	double gamma = values-alpha;
#if 0
	double chi = 1/(sigma*sqrt(T))*(std::log(K/s[0])-(mu-sigma*sigma/2)*T);
	double expected = s[0]*exp(mu*T)*gsl_cdf_gaussian_P(sigma*sqrt(T)-chi,sigma)-K*gsl_cdf_gaussian_P(-chi,sigma);
	printf("expected, chi: %f, %f\n",expected,chi);
#endif

	//get the call option value of every GBM
	for(int j=0;j<N;j++){
		s[0]=szero;
		w[0]=gsl_ran_ugaussian(r);
		for(int i=1;i<=M;i++){
			w[i]=w[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
			
			s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w[i]);
		}
		values = std::max(s[M]-K,(double)0);

		gamma = values-alpha;
		alpha+=gamma/(j+1);
		file << j+1 <<" " << alpha <<"\n";

	}


	file.close();
}


int main(){
	double szero=10;
	int N=1.0e+6;
	int K=10;
	double mu=0.1;
	double sigma=0.2;
	int T=1;
	double deltat=1;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	for(int i=1;i<=5;i++)
		runTrial(r,N,K,szero,mu,sigma,T,deltat,i);
}
