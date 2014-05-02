//a list of parameters and functions used in bubble programs
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "parameters1.h"

using namespace std;

int main()
{
cout << "choose N: " << endl;
cin >> N;
N_t = (int)(points_ratio*N);
cout << "N_t = " << N_t << endl;

for (unsigned int k=0; k<N*N_t; k++)
	{
	cout << left << "(";
	for (int l=0; l<dim; l++)
		{
		cout << intCoords(k)[l];
		if (l!=(dim-1))
			{
			cout << ",";
			}
		}
	cout << ") ";
	
	if (intCoords(k)[0]!=0 and intCoords(k)[0]!=(N_t-1))
		{
		for (unsigned int m=0; m<2*dim; m++) //over neighbours in both directions
			{
			double sign = pow(-1.0,m);
			int direc = (int)m/2;
			cout << "(";
			for (int l=0; l<dim; l++)
				{
				cout << intCoords(neigh(k,direc,sign))[l];
				if (l!=(dim-1))
					{
					cout << ",";
					}
				}
			cout << ") ";
			}
		}
	cout << endl;
	}

return 0;
}
