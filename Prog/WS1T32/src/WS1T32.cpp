#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>

int main(){
	//#values
	int n=1000000;

	//interval bounds
	int a=-3,b=3;

	//initialize rng
	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);

	//file to write to
	std::ofstream file;
	file.open("task2.dat");

	//Compute unit Gaussian distribution (sigma = 1).
	//maxpx = 1/sqrt(2*pi).
	double maxpx=gsl_ran_ugaussian_pdf(0);

	for(int i=0;i<n;i++){
	//draw x,y according to rejection sampling
	double x=(b-a)*gsl_rng_uniform(r)+a;
	double y=maxpx*gsl_rng_uniform(r);

	//step 4 of algorithm
	if(y<=gsl_ran_ugaussian_pdf(x))
		file << x <<"\n";
	else
		i--;
	}

	file.close();
	gsl_rng_free(r);

	return 0;
}
