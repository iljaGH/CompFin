#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <fstream>
#define _USE_MATH_DEFINES

int f(double x){
	return 1;
}

void trapezoidal(int l){
	int n=pow(2,l)-1;
	double nodes[n];
	double weights[n];
	double integral=0;

	for(int i=0;i<n;i++){
		nodes[i]=(i+1)/(double)(n+1);
		if(i==0 || i==n-1)
			weights[i]=(double)3/2;
		else
			weights[i]=1;

		integral+=weights[i]*f(nodes[i]);
	}
	integral=integral/(n+1);

	printf("approximate integral: %f\n",integral);
}

int main(){
	int level=4;
	int n=pow(2,level)-1;
	double nodes[n];
	std::ofstream file;

	file.open("gausslegendre.dat");

	gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(n);
	for(int i=0;i<n;i++){
		double y;
		double x;
		gsl_integration_glfixed_point(0,1,i,&x,&y,t);
		for(int j=0;j<n;j++){
			double z;
			gsl_integration_glfixed_point(0,1,j,&z,&y,t);
			file << x << " " << z <<"\n";
		}
	}

	file.close();
	file.open("cc.dat");
	for(int i=0;i<n;i++){
		double x=0.5*(1-std::cos(M_PI*(i+1)/(n+1)));
		for(int j=0;j<n;j++){
			file << x << " " <<0.5*(1-std::cos(M_PI*(j+1)/(n+1)))<<"\n";
		}		
	}

	file.close();
	file.open("trapezoidal.dat");

	for(int i=0;i<n;i++){
		double x=(i+1)/(double)(n+1);
		for(int j=0;j<n;j++){
			file << x << " " <<(j+1)/(double)(n+1)<<"\n";
		}		
	}
	file.close();

	return 1;
}
