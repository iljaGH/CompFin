#include "includes.hpp"

void runTrial(gsl_rng* r,int N, int K, int szero, double mu, double sigma,
		int T, double deltat, std::ofstream& file){

	int M=T/deltat;
	double w[M+1],s[M+1];

	double Phi = gsl_cdf_gaussian_P(,sigma);
	double closedForm = s[0]*exp(mu*T)*Phi - K*phi;

	double values = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	double gamma = 0.0;

	//get the call option value of every GBM.
	for(int j=0;j<N;j++){
		s[0]=szero;
		w[0]=gsl_ran_ugaussian(r);
		for(int i=1;i<=M;i++){
			w[i]=w[i-1]+sqrt(i*deltat-(i-1)*deltat)*gsl_ran_ugaussian(r);
			
			s[i]=s[0]*exp((mu-0.5*sigma*sigma)*i*deltat+sigma*w[i]);
		}
		values=std::max(s[M]-K,0.0);
		//mean+variance
		gamma=values-alpha;
		alpha+=gamma/(j+1);
		beta=beta+gamma*gamma*j/(j+1);
	}

	double esigma=sqrt(beta/(N-1));

	file<<alpha<<" "<<esigma<<" ";

}

int main(){
	long N=1000;
	int K=10;
	double mu=0.1;
	double sigma=0.2;
	double deltaT[] = {0.2,0.4,0.5,1.0,2.0};
	int T=2;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	std::ofstream file;

	file.precision(3);
	file.setf(std::ios::scientific);

	file.open("task2.dat");

	for (N = 1000; N <= 1.0e+7; N*=2.0e+0)
	{
		file<<N<<" ";
		for (int j = 0; j<5; ++j)
		{
			runTrial(r,N,K,10,mu,sigma,T,deltaT[j],file);
		}
		file<<std::endl;
	}

	file.close();
}
