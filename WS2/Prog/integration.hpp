#ifndef INTEGRATION_HPP_
#define INTEGRATION_HPP_
#include "includes.hpp"

#if 0
double f(double x){
	return pow(x,2);
}
#endif

void montecarlo(int l,
		std::vector<double> &nodes,
		std::vector<double> &weights){
	int n=pow(2,l)-1;

	gsl_rng* r;
	r=gsl_rng_alloc(gsl_rng_mt19937);

	weights.resize(n,1.0/(n+1));
	nodes.resize(n,1.0);

	for(int i=0;i<n;i++)
		{
			nodes[i]=gsl_rng_uniform(r);
		}
}

void trapezoidal(int l,
		std::vector<double> &nodes,
		std::vector<double> &weights){
	int n=pow(2,l)-1;

	weights.resize(n,1.0/(n+1));
	weights.front() *= 3.0/2.0;
	weights.back() *= 3.0/2.0;

	nodes.resize(n,0.0);
	for (int i = 1; i <= n; ++i)
		nodes[i-1] = (double)i/(n+1);

#if 0 //TEST
	double integral=0;

	for(int i=0;i<n;++i){
		integral += weights[i]*f(nodes[i]);
	}
	integral=integral/(n+1);

	printf("approximate integral: %f\n",integral);
#endif
}

void gaussian(int l,
		std::vector<double> &nodes,
		std::vector<double> &weights){
	int n=pow(2,l)-1;

	gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(n);//intervals

	nodes.resize(n,0.0);
	weights.resize(n,0.0);

	for (int i = 0; i < n; ++i)
	{
		 gsl_integration_glfixed_point (0, 1, i, &(*(nodes.begin() + i)),
				 &(*(weights.begin() + i)), t);
	}

	gsl_integration_glfixed_table_free(t);

}

//clenshaw-curtis quadrature
void cc(int l,
		std::vector<double> &nodes,
		std::vector<double> &weights){
	int n=pow(2,l)-1;

#if 0
	double integral=0;
#endif

	weights.resize(n,0.0);
	nodes.resize(n,0.0);

	for(int i=0;i<n;i++){
		nodes[i]=0.5*(1-std::cos(M_PI*(i+1)/(n+1)));

		double sum=0;
		for(int j=1;j<=(n+1)/2;j++){
			sum+=(double)1/(2*j-1)*std::sin((2*j-1)*(i+1)*M_PI/(n+1));
		}

		weights[i]=(double)2/(n+1)*sin((i+1)*M_PI/(n+1))*sum;
#if 0 //TEST
		integral+=weights[i]*f(nodes[i]);
#endif
	}

#if 0
	printf("%f\n",integral);
#endif

}

typedef void (*Methode) (int,
		std::vector<double> &,
		std::vector<double> &);

typedef double (*Fcn) (double, void*);

double quadrature(int l,Methode nodesWeights,Fcn f, void* parameters)
{
	std::vector<double> nodes, weights;
	double result = 0.0;

	nodesWeights(l,nodes,weights);

	for(size_t i = 0; i < nodes.size(); ++i)
		result += f(nodes[i],parameters)*weights[i];

	return result;

}

#endif /*INTEGRATION_HPP_*/
