#include "includes.hpp"
#include "integration.hpp"

double f(double x){
	return x;
}

int main(){
	int level=10;
	std::vector<double> nodes, weights;
	trapezoidal(level,nodes,weights);
}
