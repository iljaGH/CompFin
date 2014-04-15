#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_rng.h>

int main(){
	//
	printf("%lf\n", (double)rand()/RAND_MAX);

	gsl_rng *r = NULL;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	printf("%lf\n", gsl_rng_uniform(r));
	gsl_rng_free(r);

	return 0;
}
