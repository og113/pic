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
typedef complex <double> comp;
typedef Eigen::SparseMatrix<comp> spmat;
typedef Eigen::VectorXcd eigvec;
typedef vector <comp> stdvec;
typedef vector <lint> intvec;
#define pi 3.14159265359
comp i(0.0,1.0);

int main()
{
complex<double> a(0.0,0.0),b,i(0.0,1.0);
cout << "Enter b: ";
cin >> b;
vector <complex<double>> u(2);
u[0] = a;
u[1] = b;

eigvec v(2);

v(0) = a + i*b;
v(1) = b;

eigvec w = Eigen::VectorXcd::Random(2);
eigvec W = Eigen::VectorXcd::Random(2);
w *= 10.0;
W *= 10.0;

cout << "w = [" << w(0) << "," << w(1) << "]" << endl;
cout << "W = [" << W(0) << "," << W(1) << "]" << endl;

cout << "v = " << "[" << v(0) << "," << v(1) << "]" << endl;

spmat m(3,3);
int non_zeros = 3;
m.reserve(non_zeros);

m.insert(0,0) = 5;
m.insert(1,1) = i*b+pow(3.0*a+2.0*(0,1)*a*b,2);
m.insert(2,2) = {9,10};
//m.insert(2,2) = {11,12};

cout << "m = " << endl;
for (unsigned long int i=0; i<3; i++)
				{
				for (unsigned long int j=0; j<3; j++)
					{
					cout << setw(7)  << setprecision(0) << fixed << m.coeffRef(i,j);
					}
				cout << endl;
				}
return 0;
}
