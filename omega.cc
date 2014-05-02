#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_poly.h>
#include "parameters4.h"

using namespace std;

typedef unsigned long int lint;
typedef complex <double> comp;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXcd eigvec;
typedef vector <comp> stdvec;
typedef vector <lint> intvec;
#define pi 3.14159265359

int main()
{
mat h(N,N);
for (lint i=0;i<N;i++)
	{
	lint x_int = intCoords(i)[1];
	vector <double> x_i(5);
	vector <double> dx_i(3);
	vector <double> Dx_i(3);
	for (int j=0; j<5; j++)
		{
		x_i[j] = real(coord( x_int-2+j ,1));
		}
	for (int j=0; j<3; j++)
		{
		dx_i[j] = x_i[j+2]-x_i[j+1];
		Dx_i[j] = (x_i[j+2]-x_i[j])/2;
		}
	h(i,i) = pow(mass,2) + 1.0/(Dx_i[1]*dx_i[0]) + 1.0/(Dx_i[1]*dx_i[1]);
	if (i!=0)
		{
		h(i,i+1) = -1.0/(dx_i[2]*pow(Dx_i[1]*Dx_i[2],0.5));	
		}
	if (i!=(N-1))
		{
		h(i,i-1) = -1.0/(dx_i[1]*pow(Dx_i[1]*Dx_i[0],0.5));		
		}
	}

using namespace Eigen;
ComplexEigenSolver<MatrixXd> eigensolver(h);
if (eigensolver.info() != Success) abort();
cout << "The eigenvalues of h are:\n" << eigensolver.eigenvalues() << endl;
cout << "Here's a matrix whose columns are eigenvectors of A \n";
//<< "corresponding to these eigenvalues:\n"
//<< eigensolver.eigenvectors() << endl;

return 0;
}
