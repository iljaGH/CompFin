#include <cstdio>
#include <cmath>

double f(double x){
	return x;
}

void trapezoidal(int l){
	int n=pow(2,l)-1;
	double nodes[n];
	double weights[n];
	double integral=0;

	for(int i=0;i<n;i++){
		nodes[i]=(i+1)/(double)(n+1);
		if(i==0 || i==n-1)
			weights[i]=(double)3/2;
		else
			weights[i]=1;

		printf("%f:%f ",weights[i], nodes[i]);

		integral+=weights[i]*f(nodes[i]);
	}
	printf("\n");
	integral=integral/(n+1);

	printf("approximate integral: %f\n",integral);
}

int main(){
	int level=1;
	trapezoidal(level);
}