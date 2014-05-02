#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <octave/oct.h>

int factorial(int f_input)
	{
	int f_result = 1;
	for (int l=0; l<f_input; l++)
		{
		f_result *= (l+1);
		}
	return f_result;
	}

using namespace std;

int main()
{
string file;
cout << "name input file in ./data." << endl;
cin >> file;
string loc_file = "./data/" + file;
double temp_X, temp_action, temp_twaction;
vector <double> X, action, twaction;
double dross;
ifstream f;
f.open((loc_file).c_str(), ios::in);
while(!f.eof()) // reads file to end of *file*, not line
	{ 
		f >> dross >> temp_X >> dross >> dross >>  temp_action >> temp_twaction >> dross;
		X.push_back(temp_X);
		action.push_back(temp_action);
		twaction.push_back(temp_twaction);
	} 
f.close();

double result = 0.0;
double alpha = X[1]-X[0];
int N = X.size()-2;
for (int i=0; i<=N; i++)
	{
	result += action[i]*pow(-1.0,N-i)*pow(X[0]+i*alpha,N)/factorial(N-i)/factorial(i)/pow(alpha,N);
	cout << "action (X=" << X[i] << ") = " << action[i] << endl;
	}
cout << "richardson result = " << result << endl;
cout << "thin wall result = " << twaction[0] << endl;

return 0;
}
