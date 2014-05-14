#include "includes.hpp"
#include "integration.hpp"
#define _USE_MATH_DEFINES


//define integrand here
double f(double x){
	return pow(x,4);
}

int main(){
	int level=4;
	std::vector<double> nodes, weights;
	cc(level,nodes,weights);
}
