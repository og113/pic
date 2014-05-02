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
unsigned int loc;
char answer = 'y';

cout << "choose N: " << endl;
cin >> N;
N_t = (int)(points_ratio*N);
cout << "N_t = " << N_t << endl;

do{
	cout << "location: " << endl;
	cin >> loc;

	cout << "coordinates are: " << endl;
	for (unsigned int k=0;k<dim;k++)
		{
		cout << intCoords(loc)[k] << " ";
		}
	cout << endl;

	cout << "again? (y/n)" << endl;
	cin >> answer;
	} while (answer == 'y');
	
cout << "print whole lattice?" << endl;
cin >> answer;
if (answer == 'y')
	{
	cout << left;
	for (unsigned int j=0; j<N; j++)
		{
		for (unsigned int k=0; k<N_t; k++)
			{
			cout << "(";
			for (unsigned int l=0;l<dim;l++)
				{
				cout << intCoords(k+j*N_t)[l] << ",";
				}
			cout << ") ";
			}
		cout << endl;
		}
	}

return 0;
}
