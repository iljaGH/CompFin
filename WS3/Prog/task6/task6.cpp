#include <cstdio>
#include <cstdlib>
#include <cmath>

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
		while(z*pow(10,epsilon)<v*pow(10,epsilon)+1){ //multiplied both by 10e12
			v=v/p;
		}
		vdc[i]=vdc[i-1]+(p+1)*v-1;
	}
}

void halton(double** halton, int n, int d){
	int prime[203];
	primes(prime,203);
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
	int n=10;
	int d=3;

	double *halt[n];
	for(int i=0;i<n;i++)
		halt[i]=new double[d];

	halton(halt,n,d);

	print2darr(halt,n,d);
		

	return 1;
}