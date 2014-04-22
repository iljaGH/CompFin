#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <cmath>

void runTrial(double mu, double sigma1, gsl_rng* r, int n, int nfile){
	std::ofstream file;
	if(nfile==1)
		file.open("task9_1.dat");
	if(nfile==2)
		file.open("task9_2.dat");
	if(nfile==3)
		file.open("task9_3.dat");
	file << "mu=" <<mu<<" sigma="<<sigma1<<"\n";

	//alpha=mean
	double alpha=mu+sigma1*gsl_ran_ugaussian(r);
	double beta=0;
	file << "1 " << " " << alpha << " 1" << "\n";
	for(int i=1;i<n;i++){
		double gamma=mu+sigma1*gsl_ran_ugaussian(r)-alpha;
		alpha+=gamma/(i+1);
		beta=beta+gamma*gamma*i/(i+1);
		double sigma=sqrt(beta/i);


		file << i+1 << " " << std::abs(mu-alpha) << " "<< std::abs(sigma1-sigma) <<"\n";

	}

	file.close();
}

int main(){

	const double PI=2*acos(0.0);
	//#values
	int n=1000000;

	//initialize rng
	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);

	//what range?
	double mu=gsl_rng_uniform(r)*10;
	double sigma1=gsl_rng_uniform(r)*10;
	double sigma2=gsl_rng_uniform(r)*10;
	double sigma3=gsl_rng_uniform(r)*10;

	runTrial(mu,sigma1,r, n, 1);
	runTrial(mu,sigma2,r, n, 2);
	runTrial(mu,sigma3,r, n, 3);
	return 0;
}
