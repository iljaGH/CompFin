#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <fstream>

int isprime(int n){
	for(int i=2;i<=sqrt(n);i++)
		if(n%i==0)
			return 0;

	return 1;
}

void primes(int primes[], int n){
	int k=0,p=0;
	while(k<n){
		p++;
		if(isprime(p)){
			primes[k]=p;
			k++;
		}
	}
}

void vandercorput(double vdc[],int n, int p, int epsilon=12){
	vdc[0]=0;
	for(int i=1;i<n;i++){
		double z=1-vdc[i-1];
		double v=1./p;
		while(z*pow(10,epsilon)<v*pow(10,epsilon)+1){
			v=v/p;
		}
		vdc[i]=vdc[i-1]+(p+1)*v-1;
	}
}

void halton(double** halton, int n, int d){
	int prime[d+3];
	primes(prime,d+3);
	double vdc[n];

	for(int i=0;i<d;i++)
	{
		vandercorput(vdc,n,prime[i+1]);
		for(int j=0;j<n;j++)
			halton[j][i]=vdc[j];
	}
}

void print2darr(double** arr, int n, int m){
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
			printf("%f ",arr[i][j]);

		printf("\n");
	}
}

int main(){

	std::ofstream file;
	file.open("halton.dat");
	file << "#x y of first 100 Halton sequence members\n";


	int n=100;
	int d=2;

	double *halt[n];
	for(int i=0;i<n;i++)
		halt[i]=new double[d];

	halton(halt,n,d);

	for(int i=0;i<n;i++)
	{
		for(int j=0;j<d;j++)
			file <<halt[i][j] << " ";

		file << "\n";
	}

	file.close();
	file.open("uniform.dat");
	file << "#x y of 100 uniform numbers\n";

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,time(NULL));

	for(int i=0;i<n;i++)
	{
			file << gsl_rng_uniform(r) << " " << gsl_rng_uniform(r)<<"\n";
	}


		

	return 1;
}