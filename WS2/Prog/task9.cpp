#include "includes.hpp"
#include "integration.hpp"
#define _USE_MATH_DEFINES

double f(double x, void * params){
	return 1.0+exp(0.5*x);
}

int main(){
	int maxlevel=16;
	double expected=2.0*exp(0.5)-1.0;
	double parameters = 0.0;

	std::ofstream file;
	file.precision(3);
	file.setf(std::ios::scientific);

	file.open("task9.dat");
	file << "#level monte-carlo trapezoidal clenshaw-curtis gauss-legendre\n";

	for(int l=1;l<=maxlevel;++l)
	{
		printf("%i\n",l);
		file << pow(2,l)-1 << " "
				<< std::abs(expected-quadrature(l,montecarlo,f,(void*)(&parameters)))/expected<<" "
				<< std::abs(expected-quadrature(l,trapezoidal,f,(void*)(&parameters)))/expected<<" "
				<< std::abs(expected-quadrature(l,cc,f,(void*)(&parameters)))/expected<<" "
				<< std::abs(expected-quadrature(l,gaussian,f,(void*)(&parameters)))/expected<<"\n";
	}

	file.close();
}
