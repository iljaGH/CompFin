#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#if 0
int isprime(int n){
	for(int i=2;i<=sqrt(n);i++)
		if(n%i==0)
			return 0;

	return 1;
}
#endif

void primes(std::vector<int> &primes, int n){
	int p=1;
	bool prime = true;

	while(primes.size() < (size_t)n){
		p++;
		prime = true;
		for(size_t i=0;i<=primes.size();i++)
		{
			if(p%primes[i] == 0)
				{
					prime = false;
					break;
				}
		}

		if(prime)
			primes.push_back(p);
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
	std::vector<int> prime;
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
