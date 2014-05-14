#include "includes.hpp"
#include "integration.hpp"


//define integrand here
double f (double x, void * params) {
  double alpha = *(double *) params;
  if(alpha){}
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

	//integrate (gauss-legendre) f over (0,1) using interval table t
	double result=gsl_integration_glfixed(&F,0,1,t); 

	double resultIter = 0.0;
	std::vector<double> nodes, weights;
	gaussian(l,nodes,weights);

	for (size_t i = 0; i < nodes.size(); ++i)
		resultIter += f(nodes[i],&alpha)*weights[i];

	printf ("result          = % .18f\n", result);
	printf ("resultIter 	 = % .18f\n", resultIter);

	gsl_integration_glfixed_table_free(t);
}
