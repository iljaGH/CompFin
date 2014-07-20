#include "includes.hpp"
#include "integration.hpp"
#define _USE_MATH_DEFINES
#define E0 2.50662823884
#define E1 -18.61500062529
#define E2 41.39119773534
#define E3 -25.44106049637
#define F0 -8.47351093090
#define F1 23.08336743743
#define F2 -21.06224101826
#define F3 3.13082909833
#define G0 0.3374754822726147
#define G1 0.9761690190917186
#define G2 0.1607979714918209
#define G3 0.0276438810333863
#define G4 0.0038405729373609
#define G5 0.0003951896511919
#define G6 0.0000321767881768
#define G7 0.0000002888167364
#define G8 0.0000003960315187

struct my_f_params {int szero; int K; int T; double mu; double sigma;};

double NormalInverseCDF(double x){
	double p=x-0.5;
	double r;
	if (fabs(p)<0.42){
		r=pow(p,2);
		return p*(((E3*r+E2)*r+E1)*r+E0)/((((F3*r+F2)*r+F1)*r+F0)*r+1);
	} else {
		if(p<0)
			r=x;
		else
			r=1-x;

		r=log(-log(r));
		r=G0+r*(G1+r*(G2+r*(G3+r*(G4+r*(G5+r*(G6+r*(G7+r*G8)))))));

		if (p<0)
			return -r;
		else return r;
	}
}

double f(double x, void * p){
	struct my_f_params * params = (struct my_f_params *)p;
    int szero = (params->szero);
    int K=(params->K);
    int T=(params->T);
    double mu=(params->mu);
    double sigma=(params->sigma);

	return std::max(szero*exp((mu-0.5*sigma*sigma)*T+sigma*sqrt(T)*NormalInverseCDF(x))-K,(double)0);
}

int main(){
	std::ofstream file;
	file.precision(3);
	file.setf(std::ios::scientific);

	int maxlevel=16;
	//K=0
	file.open("task10_0.dat");
	file << "K=0 #t monte-carlo trapezoidal clenshaw-curtis gauss-legendre\n";

	struct my_f_params params = {10,0,2,0.1,0.2};
	double expected=12.21402758160169833921071994639;

	for(int l=1;l<=maxlevel;++l)
	{
		printf("%i\n",l);
		file << pow(2,l)-1 << " "
				<< std::abs(expected-quadrature(l,montecarlo,f,(void*)(&params)))/expected<<" "
				<< std::abs(expected-quadrature(l,trapezoidal,f,(void*)(&params)))/expected<<" "
				<< std::abs(expected-quadrature(l,cc,f,(void*)(&params)))/expected<<" "
				<< std::abs(expected-quadrature(l,gaussian,f,(void*)(&params)))/expected<<"\n";
	}

	file.close();

	//K=10
	struct my_f_params params2 = {10,10,2,0.1,0.2};
	file.open("task10_10.dat");
	file << "K=10 #t monte-carlo trapezoidal clenshaw-curtis gauss-legendre\n";

	expected=2.652809490034880944992448389588518;

	for(int l=1;l<=maxlevel;++l)
		{
			file << pow(2,l)-1 << " "
					<< std::abs(expected-quadrature(l,montecarlo,f,(void*)(&params2)))/expected<<" "
					<< std::abs(expected-quadrature(l,trapezoidal,f,(void*)(&params2)))/expected<<" "
					<< std::abs(expected-quadrature(l,cc,f,(void*)(&params2)))/expected<<" "
					<< std::abs(expected-quadrature(l,gaussian,f,(void*)(&params2)))/expected<<"\n";
		}

	file.close();
}
