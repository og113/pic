//a list of parameters and functions used in bubble programs
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

typedef unsigned long int lint;
typedef complex <double> comp;
typedef Eigen::SparseMatrix<comp> spMat;
typedef Eigen::VectorXd eigVec;
typedef vector <comp> stdVec;
typedef vector <unsigned int> intVec;
#define pi 3.14159265359
comp i(0.0,1.0);

//parameters
#define dim  2 //number of euclidean dimensions
#define dimd 2.0 //number of euclidean dimension, with .0
unsigned int N = 64; //number of points in each spatial dimension
double points_ratio = 1.2; //N_t/N
unsigned int N_t = (int)points_ratio*N; //number of points in time direction
double R = 100.0; //size of bubble
double mass = 1.0; 
double epsilon = pow(10.0,-4.0); //energy difference
double L = 4*R;//total spatial size of lattice
double L_t = 2*R; //total time size of lattice
lint kappa = (int)0.0*N_t;//first kappa time steps are in minkowskian evolution

//derived quantities
double T_M = kappa*L_t/(N_t-1);//length of time in minkowkian evolution
double T_E = (1.0-kappa/(N_t-1))*L_t;//length of time in euclidean evolution
double X = mass*R; //the large thin-wall parameter, the ratio of the size of the bubble to the size of the wall
double lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0; //quartic coupling, must have dim > 1
double v =  mass*pow(lambda,-0.5); //vacuum phi
double a = L/(N-1.0); //step sizes in each spatial dimension
double b = L_t/(N_t-1.0); //step sizes in time
double rho_step = pow((dim-1)*pow(a,2.0)+pow(b,2.0),0.5); //diagonal step size
lint Totdim = pow(N,dim-1)*N_t;

//determining number of runs
double closeness = pow(10,-3.0);


//functions

//factorial function
int factorial(const int& f_input)
	{
	int f_result = 1;
	for (int l=0; l<f_input; l++)
		{
		f_result *= (l+1);
		}
	return f_result;
	}


//gives absolute value of a number
double absolute (const double& amplitude)
	{
	double abs_amplitude;
	if (amplitude > 0)
		{
		abs_amplitude = amplitude;
		}
	else
		{
		abs_amplitude = -amplitude;
		}
	return abs_amplitude;
	}


//solid angle in any dimension
double solidAngle(const int& dimension)
	{
	double solid_angle;
	if (dimension!=1 and dimension%2==0)
		{
		solid_angle = 2*pow(pi,dimension/2.0)/factorial(dimension/2 - 1);
		}
	else if (dimension%2==1)
		{
		solid_angle = pow(2.0,dimension)*pow(pi,(dimension-1.0)/2.0)*factorial((dimension-1)/2)/factorial(dimension-1);
		}
	else
		{
		cout << "solid_angle error";
		}
	return solid_angle;
	}

//gives integer (t,x,y,...) coordinates from the single location number i
intVec intCoords(const lint& locNum)
	{
	intVec intCoordVector(dim);
	intVec param(dim);
	param[dim-1] = locNum;
	intCoordVector[dim-1] = (int)locNum/(pow(N,dim-2)*N_t);
	for (int k=1; k<dim; k++)
		{
		param[dim-1-k] = param[dim-k] - intCoordVector[dim-k]*pow(N,dim-k-1)*N_t;
		if (k!=dim-1)
			{
			intCoordVector[dim-1-k] = (int)param[dim-k-1]/(pow(N,dim-k-2)*N_t);
			}
		else if (k==(dim-1))
			{
			intCoordVector[dim-1-k] = (int)param[dim-k-1];
			}
		}
	return intCoordVector;
	}

//gives values of coordinates
stdVec coords(const lint& locNum)
	{
	stdVec coordVector (dim);
	if (intCoords(locNum)[0] < kappa)
		{
		coordVector[0] = b*intCoords(locNum)[0] + i*T_E;
		}
	else
		{
		coordVector[0] = T_M + (T_E-b*(intCoords(locNum)[0]-kappa));
		}
	for (int k=1; k<dim; k++)
		{
		coordVector[k] = -L/2.0 + a*intCoords(locNum)[k];
		}
	return coordVector;
	}
	
	//another coordinate function
comp coord(const lint& locNum, const int& direction)
    {
    comp coordinate;
    if (direction == 0)
        {
        if (locNum < kappa)
            {
            coordinate = {b*intCoords(locNum)[0],T_E}; //this isn't working at the moment
            }
        else
            {
            coordinate = {T_M,T_E-b*(intCoords(locNum)[0]-kappa)};
            }
        }
    else
        {
        coordinate = -L/2.0 + a*intCoords(locNum)[direction];
        }
    return coordinate;
    }
	

long int neigh(const lint& locNum, const unsigned int& direction, const signed int& sign) //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
	{
	long int neighLocation = -1; //this is the result if there are no neighbours for the given values of the argument
	if (direction==0)
		{
		if (sign==1 and intCoords(locNum)[0]!=(N_t-1))
			{
			neighLocation = locNum+1;
			}
		if (sign==-1 and intCoords(locNum)[0]!=0)
			{
			neighLocation = locNum-1;
			}
		}
	else if (intCoords(locNum)[direction]==0 and sign==-1)
		{
		neighLocation = locNum+(N-1)*pow(N,direction-1)*N_t;
		}
	else if (intCoords(locNum)[direction]==(N-1) and sign==1)
		{
		neighLocation = locNum-(N-1)*pow(N,direction-1)*N_t;
		}
	else
		{
		neighLocation = locNum+sign*pow(N,direction-1)*N_t;
		}
	return neighLocation;
	}

//need to sort this out now intCoords is a struct.
//a test for whether or not a location is at a boundary or not. returns zero if not at a boundary, otherwise nonzero.
int edge_test(const int& locNum, const int& degree) //degree=1 only gives positive boundaries, degree=2 gives both positive and negative
	{
	int test_result = 0;
	for (int l=0; l<dim; l++)
		{
		if (degree == 1)
			{
			if (l==0 and intCoords(locNum)[l]==(N_t-1))
				{
				test_result++;
				}
			else if(l>0 and intCoords(locNum)[l]==(N-1))
				{
				test_result++;
				}
			}
		else if (degree == 2)
			{
			if (l==0 and (intCoords(locNum)[l]==(N_t-1) or intCoords(locNum)[l]==0))
				{
				test_result++;
				}
			else if(l>0 and (intCoords(locNum)[l]==(N-1) or intCoords(locNum)[l]==0))
				{
				test_result++;
				}
			}
		else
			{
			cout << "edge test failure" << endl;
			}
		}
	return test_result;
	}

//gives dt (or dx etc.) for sites
comp Dt(const lint& locNum)
	{
	comp DT;
	if (intCoords(locNum)[0]==0)
		{
		DT = (coord(locNum+1,0)-coord(locNum,0))/2.0;
		}
	if (intCoords(locNum)[0]==(N_t-1))
		{
		DT = (coord(locNum,0)-coord(locNum-1,0))/2.0;
		}
	else
		{
		DT = (coord(locNum+1,0)-coord(locNum-1,0))/2.0 ;
		}
	return DT;
	}

//gives Dt (or Dx etc.) for links
comp dt(const lint& locNum)
	{
	comp dT;
	if (intCoords(locNum)[0]<(N_t-1))
		{
		dT = coord(locNum + 1,0)-coord(locNum,0);
		}
	else
		{
		dT = 0;
		cout << "should not be calling link dt at t=(N_t-1)" << endl;
		}
	return dT;
	}


//prints a vector of strings to the terminal
void cout_string (const vector<string>& labels)
	{
	cout << left;
	for (unsigned l=0; l<labels.size(); l++)
		{
		cout << setw(12) << labels[l];
		}
	cout << endl;
	}


//asks user questions about how what they want the program to do
void askQuestions (char * loopResponse, string * parameterChoice, int * minValue, int * maxValue, int * totalLoops, char * printChoice, int * printRun, int * minRuns)
	{
	cout << "loop through a parameter? (y/n)" << endl;
	cin >> *loopResponse;
	if (*loopResponse=='y')
		{
		cout << "which parameter?: N, R, mass, epsilon, L, X." << endl;
		cin >> *parameterChoice;
		cout << "choose min and max values (respectively)." << endl;
		cin >> *minValue >> *maxValue;
		cout << "choose number of loops." << endl;
		cin >> *totalLoops;
		}
	else
		{
		cout << "print DDS matrix or -DS vector or the action (earlier) or phi or none? (m/v/a/p/n)" << endl;
		cin >> *printChoice;
		if (*printChoice=='m' or *printChoice=='v' or *printChoice=='a' or *printChoice=='p')
			{
			cout << "choose run to print" << endl;
			cin >> *printRun;
			}
		cout << "choose minimum runs" << endl;
		cin >> *minRuns;
		cout << endl;
		}
	}


//prints action early if user chooses
void earlyAction ( const char& printChoice, const int& runsCount, const int& printRun, const double& kinetic, const double& potL, const double& potE, const double& Action)
	{
	if (printChoice == 'a' and runsCount==printRun)
		{
		cout << left;
		cout << "kinetic = " << kinetic << endl;
		cout << "pot_lambda = " << potL << endl;
		cout << "pot_epsilon = " << potE << endl;
		cout << "action = " << Action << endl;
		}
	}


//once user has chosen to vary a parameter, this changes that parameter as well as all those affected by it.
void change_parameters (const double& input_parameter, const string& parameter_name)
	{
	if (parameter_name.compare("N") == 0)
		{	
		N = input_parameter;
		N_t = input_parameter*points_ratio;
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2.0)+pow(b,2.0),0.5);
		Totdim = N_t*pow(N,dim-1);
		kappa = (int)0.5*N_t;
		T_M = kappa*L_t/(N_t-1);
		T_E = (1.0-kappa/(N_t-1))*L_t;
		}
	else if (parameter_name.compare("R") == 0)
		{
		R = input_parameter;
		L = 4*R;
		L_t = 2*R;
		X = mass*R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2.0)+pow(b,2.0),0.5);
		T_M = kappa*L_t/(N_t-1);
		T_E = (1.0-kappa/(N_t-1))*L_t;
		}
	else if (parameter_name.compare("mass") == 0)
		{
		mass = input_parameter;
		X = mass*R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		}
	else if (parameter_name.compare("epsilon") == 0)
		{
		epsilon = input_parameter;
		lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		}
	else if (parameter_name.compare("L") == 0)
		{
		L = input_parameter;
		a = L/(N-1.0);
		rho_step = pow((dim-1)*pow(a,2.0)+pow(b,2.0),0.5);
		}
	else if (parameter_name.compare("X") == 0) //in this case, X takes over from mass as one of the independent variables
		{
		X = input_parameter;
		mass = X/R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		}
	else
		{
		cout << "change_parameters function doesn't change this parameter: " << parameter_name << endl;
		}
	}

void convergenceQuestions (const int& runsCount, double* runsTest, const int minRuns, const double& Action, const clock_t& Clock, clock_t* Wait, char* stopWait, char* printWait, char* printChoice, int* printRun, const double& Kinetic, const double& potL, const double& potE, const double& Twaction, bool boolWait, bool* Convergence, const double& smallLoc, double* Gamma, vector<int> Sigma)
	{
 if (*runsTest<=closeness and runsCount>=minRuns) //if on last run
	{
//crude test for convergence
	if (absolute(Action) < absolute(Twaction)*pow(10,3))
		{
		*Convergence = true;
		}
	}

//test to see if action grows large and how large bubble is
	if (absolute(Action) > absolute(Twaction)*pow(10,3))
		{
		*Convergence = false;
		*runsTest = 0;
		double gamma_default = 0.01;
		unsigned int l=1;
		unsigned int n = 1 + Sigma.size();
		if (n>10)
			{
			cout << "could not get bubble to stay still" << endl;
			*Convergence = true;
			}
		else if (smallLoc<0.8*R or smallLoc>1.2*R)
			{
			if (smallLoc<0.8*R)
				{
				Sigma.push_back(1);
				}
			if (smallLoc>1.2*R)
				{
				Sigma.push_back(-1);
				}
			if (n>1)
				{
				for (unsigned int j=1; j<n;j++)
					{
					if ( Sigma[0]*Sigma[j]==1 )
						{
						l++; //l gives first time sign differs from original sign, this dictates the max size of changes in gamma
						}
					else
						{
						break;
						}
					}
				}
			else
				{
				l=1;
				}
			*Gamma += gamma_default*Sigma[n-1]*pow(2.0,2.0*l-n-1.0);
			}
		else
			{
//			cout << "action grew large and not sure what happened to bubble so stopped loop" << endl;
			*runsTest = 0;	
			}
		}

	//stopping newton-raphson loop if action doesn't converge after 1000 loops
	if (runsCount > 1000)
		{
		cout << "over 1000 runs without convergence of action to " << closeness << ", increase value of 'closeness' to get a result" << endl;
		*runsTest = 0;
		}

//quick test to make sure newton-raphson loop doesn't end early
	if (Action == 1.0)
		{
		cout << "action = 2, change initial value of action to prevent while loop from ending" << endl;
		}

//prints runs_count if looping is taking ages
	if((Clock-*Wait)/1000000.0>300 and not(boolWait))
		{
		cout << "number of newton-raphson loops = " << runsCount << endl;
		*Wait = Clock;
		cout << "stop this looping and move on? (y/n)" << endl;
		cin >> *stopWait;
		if ( *stopWait == 'y')
			{
			*runsTest = 0;
			*Convergence = true;
			}
		else
			{
			cout << "print phi and action on the next loop? (y/n)" << endl;
			cin >> *printWait;
			if (*printWait == 'y')
				{
				*printChoice = 'p';
				*printRun = runsCount+1;
				cout << left;
				cout << "kinetic = " << Kinetic << endl;
				cout << "pot_lambda = " << potL << endl;
				cout << "pot_epsilon = " << potE << endl;
				cout << "action = " << Action << endl;
				}
			}
		}
	}

