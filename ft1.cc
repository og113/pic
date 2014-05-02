#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_poly.h>
#include "parameters1.h"

using namespace std;

int main()
{
Eigen::MatrixXd h(N,N);
eigVec p(Totdim), FTp(Totdim);
int n=0; //fourier transform with w_k given by the nth eigenvalue of h

//alpha gives the range over which the thin wall solution tanh is used
double alpha = 20.0;
//gamma gives an extra bit of size to the bubble
double gamma = 0.0;

//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
double b = -pow(v,2.0);
double c = epsilon/v/lambda;
vector<double> root(3);
gsl_poly_solve_cubic (0, b, c, &root[0], &root[1], &root[2]);
sort(root.begin(),root.end());

//assigning thin wall solution as initial guess for phi
for (unsigned long int i = 0; i<Totdim; i++)
	{
	double bubble_size = R * (1.0 + gamma);
	vector <long unsigned int> int_x = coords(i);
	vector <double> x(dim);
	x[0] = -L/beta + a*int_x[0]/beta; //time coordinate steps are 1/beta the size
	double rho_sqrd = pow(x[0],2.0);
	for (int j=1; j<dim; j++)
		{
		x[j] = -L/2.0 + a*int_x[j];
		rho_sqrd += pow(x[j],2.0);
		}
	double rho = pow(rho_sqrd,0.5);
	if (R<alpha/mass)
		{
		cout << "X = R*mass is too small. It should not be less than " << alpha << endl;
		}
	else
		{
//		p(i) = root[0];
		if (rho<(bubble_size-alpha/mass))
			{
			p(i) = root[0];
			}
		else if (rho>(bubble_size+alpha/mass))
			{
			p(i) = root[2];
			}
		else
			{
			p(i) = v*tanh(mass*(rho-bubble_size)/2.0);
			}
		}
	}


h.setZero();
for (lint i=0;i<N;i++)
	{
	if (i>0 and i<(N-1))
		{
		h(i,i) = pow(mass,2) + 1.0/pow(a,2) + 1.0/pow(a,2);
		h(i,i+1) = -1.0/pow(a,2);
		h(i,i-1) = -1.0/pow(a,2);
		}
	else if (i==0)
		{
		h(i,i) = pow(mass,2) + 1.0/pow(a,2) + 1.0/pow(a,2);
		h(i,i+1) = -1.0/pow(a,2);
		}
	else
		{
		h(i,i) = pow(mass,2) + 1.0/pow(a,2) + 1.0/pow(a,2);
		h(i,i-1) = -1.0/pow(a,2);
		}
	}

//solving for eigenvalues and eigenvectors
Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(h);
if (eigensolver.info() != Eigen::Success) abort();

//pair of eigenvalue and eigenvector
struct ePair{
complex <double> eVal;
Eigen::VectorXcd eVec;
};

//clever swapper
struct by_eVal {
	bool operator()(ePair const &eigen1, ePair const &eigen2)
	{ 
    return abs(eigen1.eVal) < abs(eigen2.eVal);
    }
};

vector <ePair> h_eigs(N);

for (unsigned int i=0; i<N; i++)
	{
	h_eigs[i].eVal = eigensolver.eigenvalues()[i];
//	eig_guess[i] = pow(pow(mass,2) + pow(2.0*pi*(i+1.0)/L,2),0.5);
	h_eigs[i].eVec = eigensolver.eigenvectors().col(i);
	}

//sort eigs
sort (h_eigs.begin(),h_eigs.end(), by_eVal());

bool print = false;
if (print)
	{
	for (unsigned int i=0; i<N; i++)
		{
		cout << h_eigs[i].eVal <<endl;
		}
	}

for (lint i=0;i<N_t;i++)
	{
	for (lint j=0;j<N;j++)
		{
		FTp(i) = pow(a,0.5)* h_eigs[n].eVec(j) * p(i + N_t*j);
		}
	}

return 0;
}
