#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>

int main(){
	/*
	 * Use standard C function to create uniform distributed random
	 * variables. Force conversion to double (return type of rand()
	 * is int.
	 */
	printf("%lf\n", (double)rand()/RAND_MAX);

	/*
	 * Create uniform distributed random variable with GSL.
	 */
	gsl_rng *r = NULL;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	printf("%lf\n", gsl_rng_uniform(r));
	gsl_rng_free(r);

	return 0;
}
