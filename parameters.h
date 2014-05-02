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
typedef vector <unsigned int> intVec;
typedef vector <comp> stdVec;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd eigVec;
#define pi 3.14159265359
#define im(argument) imag(argument)
#define re(argument) real(argument)
comp i(0.0,1.0);

//parameters
#define dim  2 //number of euclidean dimensions
#define dimd 2.0 //number of euclidean dimension, with .0
unsigned int N = 64; //number of points in each spatial dimension
double points_ratio = 1.0; //N_t/N
unsigned int N_t = (int)(points_ratio*N); //number of points in euclidean time direction
unsigned int N_t_m = N_t;
double R = 100.0; //size of bubble
double mass = 1.0; 
double epsilon = pow(10.0,-4.0); //energy difference
double L_t = 1.2*R; //length of euclidean time path

//derived quantities
double theta = asin(L_t/2.0/R); //angle at which circles meet
double L = 1.5*L_t;//total spatial size of lattice - 1.5 gives leeway
double X = mass*R; //the large thin-wall parameter, the ratio of the size of the bubble to the size of the wall
double lambda = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/R/3.0; //quartic coupling, must have dim > 1
double v =  mass*pow(lambda,-0.5); //vacuum phi
double a = L/(N-1.0); //step sizes in each spatial dimension
double b = L_t/(N_t-1.0); //step sizes in time
double rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5); //diagonal step size
double L_t_m = b*N_t_m;
lint Eucdim = pow(N,dim-1)*N_t;
lint Minkdim = pow(N,dim-1)*N_t_m;
lint Totdim = pow(N,dim-1)*(N_t+N_t_m);

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
	if (dimension%2==0)
		{
		solid_angle = 2*pow(pi,dimension/2)/factorial(dimension/2 - 1);
		}
	else if (dimension%2==1)
		{
		solid_angle = pow(2.0,dimension)*pow(pi,(dimension-1)/2)*factorial((dimension-1)/2)/factorial(dimension-1);
		}
	else
		{
		cout << "solid_angle error";
		}
	return solid_angle;
	}


//gives integer (t,x,y,...) coordinates from the single location number i
intVec intCoords(const lint& locNum, const unsigned int& Nt)
	{
	intVec intCoordVector(dim);
	intVec param(dim);
	param[dim-1] = locNum;
	intCoordVector[dim-1] = (int)locNum/(pow(N,dim-2)*Nt);
	for (int k=1; k<dim; k++)
		{
		param[dim-1-k] = param[dim-k] - intCoordVector[dim-k]*pow(N,dim-k-1)*Nt;
		if (k!=dim-1)
			{
			intCoordVector[dim-1-k] = (int)param[dim-k-1]/(pow(N,dim-k-2)*Nt);
			}
		else if (k==(dim-1))
			{
			intCoordVector[dim-1-k] = (int)param[dim-k-1];
			}
		}
	return intCoordVector;
	}

//gives values of coordinates in euclidean part
eigVec coords(const lint& locNum)
	{
	eigVec coordVector (dim);
	coordVector[0] = i*(-L_t + b*intCoords(locNum,N_t)[0]);
	for (int k=1; k<dim; k++)
		{
		coordVector[k] = -L/2.0 + a*intCoords(locNum,N_t)[k];
		}
	return coordVector;
	}
	
	//another coordinate function in euclidean part
comp coord(const lint& locNum, const int& direction)
    {
    comp coordinate;
    if (direction == 0)
        {
		coordinate = i*(-L_t + b*intCoords(locNum,N_t)[0]);
        }
    else
        {
        coordinate = -L/2.0 + a*intCoords(locNum,N_t)[direction];
        }
    return coordinate;
    }
	

lint neigh(const lint& locNum, const unsigned int& direction, const signed int& sign, const unsigned int& Nt) //periodic in space but not time, note this does not give any neighbours if intCoords(locNum)[0]==0 or N_t-1
	{
	long int neighLocation = -1; //this is the result if there are no neighbours for the given values of the argument
	if (direction==0)
		{
		if (sign==1 and intCoords(locNum,Nt)[0]!=(Nt-1))
			{
			neighLocation = locNum+1;
			}
		if (sign==-1 and intCoords(locNum,Nt)[0]!=0)
			{
			neighLocation = locNum-1;
			}
		}
	else if (intCoords(locNum,Nt)[direction]==0 and sign==-1)
		{
		neighLocation = locNum+(N-1)*pow(N,direction-1)*Nt;
		}
	else if (intCoords(locNum,Nt)[direction]==(N-1) and sign==1)
		{
		neighLocation = locNum-(N-1)*pow(N,direction-1)*Nt;
		}
	else
		{
		neighLocation = locNum+sign*pow(N,direction-1)*Nt;
		}
	return neighLocation;
	}

//need to sort this out now intCoords is a struct.
//a test for whether or not a location is at a boundary or not. returns zero if not at a boundary, otherwise nonzero.
int edge_test(const int& locNum, const int& degree, const unsigned int& Nt) //degree=1 only gives positive boundaries, degree=2 gives both positive and negative
	{
	int test_result = 0;
	for (int l=0; l<dim; l++)
		{
		if (degree == 1)
			{
			if (l==0 and intCoords(locNum,Nt)[l]==(Nt-1))
				{
				test_result++;
				}
			else if(l>0 and intCoords(locNum,Nt)[l]==(N-1))
				{
				test_result++;
				}
			}
		else if (degree == 2)
			{
			if (l==0 and (intCoords(locNum,Nt)[l]==(Nt-1) or intCoords(locNum,Nt)[l]==0))
				{
				test_result++;
				}
			else if(l>0 and (intCoords(locNum,Nt)[l]==(N-1) or intCoords(locNum,Nt)[l]==0))
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
	
//gives time step for links
comp dt (const lint & locNum)
	{
	comp dtime;
	if (intCoords(locNum,N_t)[0]==(N_t-1))
		{
		dtime = 0;
		}
	else
		{
		dtime = coord(locNum+1,0) - coord(locNum,0);
		}
	return dtime;
	}

//gives time step for sites
comp Dt (const lint & locNum)
	{
	comp Dtime;
	if (intCoords(locNum,N_t)[0]==0)
		{
		Dtime = (coord(locNum+1,0) - coord(locNum,0))/2.0;
		}
	else if (intCoords(locNum,N_t)[0]==(N_t-1))
		{
		Dtime = (coord(locNum,0) - coord(locNum-1,0))/2.0;
		}
	else
		{
		Dtime = (coord(locNum+1,0) - coord(locNum-1,0))/2.0;
		}
	return Dtime;
	}
	

//prints a vector of strings to the terminal
void coutStringShort (const vector<string>& labels)
	{
	cout << left;
	for (unsigned l=0; l<labels.size(); l++)
		{
		cout << setw(8) << labels[l];
		}
	}
	
//prints a vector of strings to the terminal
void coutStringLong (const vector<string>& labels)
	{
	cout << left;
	for (unsigned l=0; l<labels.size(); l++)
		{
		cout << setw(20) << labels[l];
		}
	cout << endl;
	}


//asks user questions about how what they want the program to do
void askQuestions (char * inputP, char * perturbResponse, char * loopResponse, string * parameterChoice, int * minValue, int * maxValue, int * totalLoops, char * printChoice, int * printRun, int * minRuns)
	{
	cout << "spherical bubble, periodic instanton, true vacuum or false vacuum (b,p,t,f)?" << endl;
	cin >> *inputP;
	cout << "add small (10^-4) perturbations in real, imaginary or both parts of the field? (r/i/b/n)" << endl;
	cin >> *perturbResponse;
	cout << "loop through a parameter? (y/n)" << endl;
	cin >> *loopResponse;
	if (*loopResponse=='y')
		{
		cout << "which parameter?: N, R, mass, epsilon, L_t, X, lambda." << endl;
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
void earlyAction ( const char& printChoice, const int& runsCount, const int& printRun, const comp& kinetic, const comp& potL, const comp& potE, const comp& Action)
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
	
void earlyVector (const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
			for (unsigned long int j=0; j<Eucdim; j++)
				{
				F << left;
				for (int r=0; r<dim; r++)
					{
					F << setw(15) << re(coord(j,r)) << setw(15) << im(coord(j,r));
					}
				F << setw(15) << vecToPrint(2*j) << setw(15) << vecToPrint(2*j+1)  << endl;		
				}
			F.close();
	}


//once user has chosen to vary a parameter, this changes that parameter as well as all those affected by it.
void changeParameters (const double& inputParameter, const string& parameterName, const char& inputP)
	{
	if (inputP == 'b')
		{
		L_t = 2*R;
		L = 4*R;
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
		
		if (parameterName.compare("N") == 0)
			{
			N = inputParameter;
			N_t = (int)(inputParameter*points_ratio);
			a = L/(N-1.0);
			b = L_t/(N_t-1.0);
			rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
			Eucdim = pow(N,dim);
			}
		else if (parameterName.compare("R") == 0)
			{
			R = inputParameter;
			L = 4*R;
			L_t = 2*R;
			X = mass*R;
			lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
			v =  mass*pow(lambda,-0.5);
			a = L/(N-1.0);
			b = L_t/(N_t-1.0);
			rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
			}
		else if (parameterName.compare("mass") == 0)
			{
			mass = inputParameter;
			X = mass*R;
			lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
			v =  mass*pow(lambda,-0.5);
			}
		else if (parameterName.compare("epsilon") == 0)
			{
			epsilon = inputParameter;
			lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
			v =  mass*pow(lambda,-0.5);
			}
		else if (parameterName.compare("L") == 0)
			{
			L = inputParameter;
			a = L/(N-1.0);
			rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
			}
		else if (parameterName.compare("L_t") == 0)
			{
			L_t = inputParameter;
			b = L_t/(N_t-1.0);
			rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
			}
		else if (parameterName.compare("X") == 0) //in this case, X takes over from mass as one of the independent variables
			{
			X = inputParameter;
			mass = X/R;
			lambda = 2.0*(dimd-1.0)*pow(mass,3.0)/epsilon/R/3.0;
			v =  mass*pow(lambda,-0.5);
			}
		}
	else if (parameterName.compare("N") == 0)
		{	
		N = inputParameter;
		N_t = (int)(inputParameter*points_ratio);
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
		Eucdim = N_t*pow(N,dim-1);
		}
	else if (parameterName.compare("R") == 0)
		{
		R = inputParameter;
		L_t = 1.2*R;
		theta = asin(L_t/2/R);
		L = 1.5*L_t*tan(theta);
		X = mass*R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
		}
	else if (parameterName.compare("mass") == 0)
		{
		mass = inputParameter;
		X = mass*R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		}
	else if (parameterName.compare("epsilon") == 0)
		{
		epsilon = inputParameter;
		lambda = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		}
	else if (parameterName.compare("L_t") == 0)
		{
		L_t = inputParameter;
		R = L_t/1.2;
		theta = asin(L_t/2/R);
		L = 1.5*L_t*tan(theta);
		X = mass*R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
		}
	else if (parameterName.compare("X") == 0) //in this case, X takes over from mass as one of the independent variables
		{
		X = inputParameter;
		mass = X/R;
		lambda = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/R/3.0;
		v =  mass*pow(lambda,-0.5);
		}
	else if (parameterName.compare("lambda") == 0) //to check 1/lambda scaling, epsilon also scales
		{
		epsilon *= lambda/inputParameter; //epsilon scales like 1/lambda
		lambda = inputParameter;
		R = 2.0*(dimd-1.0)*pow(mass,3)/epsilon/lambda/3.0;
		X = mass*R;
		v =  mass*pow(lambda,-0.5);
		L_t = 1.2*R;
		theta = asin(L_t/2.0/R);
		L = 1.5*L_t*tan(theta);
		a = L/(N-1.0);
		b = L_t/(N_t-1.0);
		rho_step = pow((dim-1)*pow(a,2)+pow(b,2),0.5);
		}
	else
		{
		cout << "change_parameters function doesn't change this parameter: " << parameterName << endl;
		}
	}

void convergenceQuestions (const int& runsCount, double* runsTest, const int minRuns, const comp& Action, const clock_t& Clock, clock_t* Wait, char* stopWait, char* printWait, char* printChoice, int* printRun, const comp& Kinetic, const comp& potL, const comp& potE, const double& Twaction, bool boolWait)
	{
	//stopping newton-raphson loop if action doesn't converge after 1000 loops
	if (runsCount > 1000)
		{
		cout << "over 1000 runs without convergence of action to " << closeness << ", increase value of 'closeness' to get a result" << endl;
		*runsTest = 0;
		}

//quick test to make sure newton-raphson loop doesn't end early
	if (Action == 2.0)
		{
		cout << "action = 2, change initial value of action to prevent while loop from ending" << endl;
		*runsTest = 0;
		}

//prints runs_count if looping is taking ages
	if((Clock-*Wait)/1000000.0>600 and not(boolWait))
		{
		cout << "number of newton-raphson loops = " << runsCount << endl;
		*Wait = Clock;
		cout << "stop this looping and move on? (y/n)" << endl;
		cin >> *stopWait;
		if ( *stopWait == 'y')
			{
			*runsTest = 0;
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


