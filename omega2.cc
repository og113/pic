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

int main()
{
Eigen::MatrixXd h(N,N);
Eigen::MatrixXcd omega(N,N);
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

vector <double> eig_guess(N);

//clever swapper
struct by_eVal {
	bool operator()(ePair const &eigen1, ePair const &eigen2) { 
    return abs(eigen1.eVal) < abs(eigen2.eVal);
    }
};

vector <ePair> h_eigs(N);

for (unsigned int i=0; i<N; i++)
	{
	h_eigs[i].eVal = eigensolver.eigenvalues()[i];
	eig_guess[i] = pow(mass,2) + pow(2.0*sin(pi*i/N)/a,2.0);
	h_eigs[i].eVec = eigensolver.eigenvectors().col(i);
	}

//sort eigs
sort (h_eigs.begin(),h_eigs.end(), by_eVal());
sort (eig_guess.begin(),eig_guess.end());

bool print = true;
if (print)
	{
	for (unsigned int i=0; i<N; i++)
		{
		cout << h_eigs[i].eVal << " " << eig_guess[i] << endl;
		}
	}

for (lint i=0;i<N;i++)
	{
	for (lint j=0;j<N;j++)
		{
		for (lint k=0;k<N;k++)
			{
			omega(i,j) = a*h_eigs[k].eVal * h_eigs[k].eVec(i) * h_eigs[k].eVec(j);
			}
		}
	}

return 0;
}
