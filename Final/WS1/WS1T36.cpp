#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <cmath>

int main(){

	const double PI=2*acos(0.0);
	//#values
	int n=1000;

	//initialize rng
	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);

	//file to write to
	std::ofstream file;
	file.open("task6.dat");

	for(int i=0;i<n;i++){
	//draw u1,u2 uniform on [0,1]
	double u1=gsl_rng_uniform(r);
	double u2=gsl_rng_uniform(r);

	double z1=sqrt(-2*log(u1))*cos(2*PI*u2);
	double z2=sqrt(-2*log(u1))*sin(2*PI*u2);

	//write to file
	file << z1<< " " << z2 <<"\n";
	}

	file.close();
	gsl_rng_free(r);

	return 0;
}
