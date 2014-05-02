#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_poly.h>

using namespace std;

typedef unsigned long int lint;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;

int main()
{
double kappa = 0.2;
double L=100.0;
unsigned int N = 1000;
double b = L/(N-1);
vec x(N), x_h(N), v(N), a(N), a_h(N);
x(0) = 1.0;
v(0) = 0.0;
a(0) = -kappa*x(0);
x_h(0) = x(0) + v(0)*b/2.0 + a(0)*pow(b,2)/4.0;
a_h(0) = -kappa*x_h(0);
v(1) = v(0) + a_h(0)*b;
x(1) = x_h(0) + v(1)*b/2.0;
for (unsigned int j=1; j<(N-1); j++)
	{
	x_h(j) = x(j) + v(j)*b/2.0;
	a_h(j) = -kappa*x_h(j);
	v(j+1) = v(j) + a_h(j)*b;
	x(j+1) = x_h(j) + v(j+1)*b/2.0;
	}

//printing output phi
	ofstream f;
	f.open("./data/lf.dat");
	for (unsigned int j=0; j<N; j++)
		{
		f << left;
		f << setw(15) << j*b ;
		f << setw(15) << x(j) << setw(15) << v(j) << endl;	
		}
	f.close();

return 0;
}
