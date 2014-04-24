#include <gsl/gsl_integration.h>
#include <cmath>
#include <cstdio>

//define integrand here
double f (double x, void * params) {
  double alpha = *(double *) params;
  double f = 2*pow(x,5)+pow(x,3)+4*pow(x,2)+3;
  return f;
}

int main ()
{
	//level
	int l=2;
  	gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(pow(2,l)-1);//intervals
  
	double alpha = 1.0; //weights

	gsl_function F;
	F.function = &f;
	F.params = &alpha;

	//integrate (gauß-legendre) f over (0,1) using interval table t
	double result=gsl_integration_glfixed(&F,0,1,t); 

	printf ("result          = % .18f\n", result);

	gsl_integration_glfixed_table_free(t);
}