#include <gsl/gsl_integration.h>
#include <cmath>
#include <cstdio>
#define _USE_MATH_DEFINES


//define integrand here
double f(double x){
	return pow(x,4);
}

//clenshaw-curtis quadrature
void cc(int l){
	int n=pow(2,l)-1;
	double nodes[n];
	double weights[n];
	double integral=0;

	for(int i=0;i<n;i++){
		nodes[i]=0.5*(1-std::cos(M_PI*(i+1)/(n+1)));

		double sum=0;
		for(int j=1;j<=(n+1)/2;j++){
			sum+=(double)1/(2*j-1)*std::sin((2*j-1)*(i+1)*M_PI/(n+1));
		}

		weights[i]=(double)2/(n+1)*sin((i+1)*M_PI/(n+1))*sum;

		integral+=weights[i]*f(nodes[i]);
	}

	printf("%f\n",integral);

}

int main(){
	int level=4;
	cc(level);
}