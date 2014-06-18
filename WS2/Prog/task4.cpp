#include "includes.hpp"

	double szero=10;
	int N=1.0e+6;
	int K=10;
	double mu=0.1;
	double sigma=0.2;
	int T=1;
	double deltat=1;


double BlackScholes(double sigma){
	double r=mu;
	double d1=(log(szero/K)+(r+sigma*sigma*0.5)*T)/(sigma*sqrt(T));
	double d2=d1-sigma*sqrt(T);
	return szero*gsl_cdf_ugaussian_P(d1)-K*exp(-r*T)*gsl_cdf_ugaussian_P(d2);
}

void runTrial(gsl_rng* r,int N, int K, int szero, double mu, double sigma, int T, double deltat, int trial){

	std::ofstream file,file2;

	switch(trial){
		case 1:file.open("task4_1.dat");file2.open("task4_1_conv.dat");break;
		case 2:file.open("task4_2.dat");file2.open("task4_2_conv.dat");break;
		case 3:file.open("task4_3.dat");file2.open("task4_3_conv.dat");break;
		case 4:file.open("task4_4.dat");file2.open("task4_4_conv.dat");break;
		case 5:file.open("task4_5.dat");file2.open("task4_5_conv.dat");break;
	}

	int M=T/deltat;
	double w[M+1],s[M+1];
	double expected=BlackScholes(sigma);

	double sum=0;
	/*double chi = 1/(sigma*sqrt(T))*(std::log((double)K/szero)-(mu-sigma*sigma/2)*T);
	double expected = szero*exp(mu*T)*gsl_cdf_gaussian_P(sigma*sqrt(T)-chi,sigma)-K*gsl_cdf_gaussian_P(-chi,sigma);
	printf("expected, chi: %f, %f\n",expected,chi);*/

	//get the call option value of every GBM
	for(int j=0;j<N;j++){
		s[0]=szero;
		w[0]=0;
		for(int i=1;i<=M;i++){
			w[i]=w[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
			
			s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w[i]);
		} 

		sum=(sum*j+std::max(s[M]-K,(double)0))/(j+1);
		file << j+1 <<" " << exp(-mu*T)*sum <<"\n";
		file2 << j+1 <<" " << std::abs(expected-exp(-mu*T)*sum)/expected<<"\n";

	}


	file.close();
	file2.close();
}

int main(){

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	for(int i=1;i<=5;i++)
		runTrial(r,N,K,szero,mu,sigma,T,deltat,i);
}
