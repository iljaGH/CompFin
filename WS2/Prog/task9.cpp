#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#define _USE_MATH_DEFINES

double f(double x, void * params){
	return 1+0.1*exp(0.5*x);
}

double cc(int l){
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

		integral+=weights[i]*f(nodes[i],0);
	}

	return integral;

}

double trapezoidal(int l){
	int n=pow(2,l)-1;
	double nodes[n];
	double weights[n];
	double integral=0;

	for(int i=0;i<n;i++){
		nodes[i]=(i+1)/(double)(n+1);
		if(i==0 || i==n-1)
			weights[i]=((double)3/2)/(n+1);
		else
			weights[i]=1./(n+1);

		printf("%f:%f ",nodes[i],weights[i]);

		integral+=weights[i]*f(nodes[i],0);
	}

	return integral;
}

double montecarlo(int l){
	int n=pow(2,l)-1;
	double sum=0;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);

	for(int i=0;i<n;i++)
		sum+=f(gsl_rng_uniform(r),0);
	
	return sum/n;
}

int main(){
	int maxlevel=10;
	double expected=2.29744;

	std::ofstream file;
	file.open("task9.dat");
	file << "#level monte-carlo trapezoidal clenshaw-curtis gauss-legendre\n";
	for(int l=2;l<=2;l++)
	{
		//Gauss-Legendre
		gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(pow(2,l)-1);//intervals
		double alpha = 1.0; //weights

		gsl_function F;
		F.function = &f;
		F.params = &alpha;

		//file << l << " "<<std::abs(expected-montecarlo(l))/expected<<" "<< std::abs(expected-trapezoidal(l))/expected<<" "<<std::abs(expected-cc(l))/expected<<" "<<std::abs(expected-gsl_integration_glfixed(&F,0,1,t))/expected<<"\n";
		gsl_integration_glfixed_table_free(t);
		printf("%f ",cc(l));
	}

	file.close();
}