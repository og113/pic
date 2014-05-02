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
typedef Eigen::SparseMatrix<comp> spMat;
typedef Eigen::VectorXcd eigVec;
typedef vector <comp> stdVec;
typedef vector <lint> intVec;
#define pi 3.14159265359

int main()
{
//clear action output file
//{
//ofstream actionfile;
//actionfile.open ("./data/actionloop.dat", ios::out | ios::trunc );
//actionfile.close();
//}

//quantities carrying user responses
char loop_response;
string parameter_choice;
int min_value;
int max_value;
int total_loops;
char print_choice;
int print_run;
int min_runs;

//asking user questions
askQuestions (&loop_response, &parameter_choice, &min_value, &max_value, &total_loops, &print_choice, &print_run, &min_runs);

//defining a time
clock_t time;
clock_t wait;

//double angle in dim
double solid_angle = solidAngle(dim);


//if user chose no loops, set total_loops to 1
if (loop_response=='n')
	{
	total_loops = 1;
	}
//begin loop over values of parameter user chose to vary
for (int loop=0; loop<total_loops; loop++)
{
//giving values of varying parameters
if (parameter_choice.compare("N") == 0 and loop_response=='y')
	{
	min_runs = 10;
	N = min_value + (int)(max_value - min_value)*loop/(total_loops-1);
	change_parameters (N, "N");
	}
else if (loop_response=='y')
	{
	min_runs = 10;
	double loop_parameter = min_value + (max_value - min_value)*loop/(total_loops-1.0);
	change_parameters (loop_parameter,parameter_choice);
	}

//starting the clock
time = clock();
wait = clock();

//defining some important scalar quantities
complex <long double> log_det_DDS;
comp action = 2.0;
double S_1 = epsilon*R/(dimd-1);
double twaction = -solid_angle*epsilon*pow(R,dimd)/dimd + solid_angle*pow(R,dimd-1.0)*S_1;
double real_small_loc = 0;
double imag_small_loc = 0;
comp eigen_zero;
//alpha gives the range over which the thin wall solution tanh is used
double alpha = 20.0;
//gamma gives an extra bit of size to the bubble
double gamma = 0.0;

//defining some quantities used to stop the Newton-Raphson loop when action stops varying
comp action_last = 1.0;
int runs_count = 0;
double runs_test = 1.0;

bool convergence = false;
do {
//initializing phi (=p)
eigVec p(Totdim);
eigVec neg(Totdim); //negative mode

#define pr(j) p(2*j)
#define pi(j) p(2*j+1)

//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
double b = -pow(v,2.0);
double c = epsilon/v/lambda;
vector<double> root(3);
gsl_poly_solve_cubic (0, b, c, &root[0], &root[1], &root[2]);
sort(root.begin(),root.end());

//assigning thin wall solution as initial guess for phi
for (lint j = 0; j<Totdim; j++)
	{
	double bubble_size = R * (1.0 + gamma);
	double rho_sqrd = 0;
	for (int k=0; k<dim; k++)
		{
		rho_sqrd += norm(coords(intCoords(j)[k]));
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
			p(j) = root[0];
			neg(j) = 0;
			}
		else if (rho>(bubble_size+alpha/mass))
			{
			p(j) = root[2];
			neg(j) = 0;
			}
		else
			{
			p(j) = v*tanh(mass*(rho-bubble_size)/2.0);
			neg(j) = -v*mass/(2*pow(cosh(mass*(rho-bubble_size)/2.0),2.0));
			}
		}
	}
	
//aught to begin another loop here whereby i loop over increasing values of T_E




//beginning newton-raphson loop, continuing loop until ratio of change of action to action is smaller than predefined "closeness".
while (runs_test>closeness or runs_count<min_runs)
{
	//quantities used to stop newton-raphson loop
	runs_test = abs(action - action_last)/abs(action_last);
	runs_count ++;
	action_last = action;

	// allocating memory for DS, DDS and delta
	int non_zeros = (Totdim*(1-2/N))*(2*dim+1)+ (2*Totdim/N)*(Totdim);//number of non-zero elements in DDS - not sure if there are too many - this assumes the time boundaries are dense
	eigVec minusDS(Totdim);
	spMat DDS(Totdim,Totdim);
	DDS.reserve(non_zeros);
	eigVec delta(Totdim);  //delta=newphi-phi in Lowdim
	eigVec dp_dx(Totdim); //derivative of phi in the x direction, a zero mode

	//setting action to zero, and defining quantities which split up action into its parts
	action = 0.0; //defined as minkowskian
	comp kinetic = 0.0;
	comp pot_l = 0.0;
	comp pot_e = 0.0;
	double real_small = real(p(0)); //is this the right thing to consider, not the abs?
	double imag_small = imag(p(0));
	//assigning values to minusDS and DDS and evaluating action
	for (unsigned long int j = 0; j < Totdim; j++)
		{
		double spatialMeasure = pow(a,dim-1);
		comp siteMeasure = spatialMeasure*Dt(j); //measure for sites in time
		stdVec linkMeasure = spatialMeasure*dt(j); //measure for links in time
		minusDS(j) = 0.0;//initializing to zero
		pot_l += -siteMeasure*lambda*pow(pow(p(j),2.0)-pow(v,2.0),2.0)/8.0;
		pot_e += -siteMeasure*epsilon*(p(2*j)-v)/v/2.0;
		if (intCoords(j)[0]!=(N_t-1))
			{
	      		for (unsigned int k=0; k<dim; k++) //just over neighbours in positive directions //need to change edge test for spatial boundaries, and better to use neighbours than this ad-hoc thing
	           		{
					if(k==0)
						{			
						kinetic += linkMeasure*pow(pr(neigh(j,k,1))+i*pi(neigh(j,k,1))-pr(j)-i*pi(j),2.0)/pow(dt(j),2.0)/2.0;
						}
					else
						{
						kinetic += -siteMeasure*pow(pr(neigh(j,k,1))+i*pi(neigh(j,k,1))-pr(j)-i*pi(j),2.0)/pow(a,2.0)/2.0;
						}
	           		}
			}
		if (intCoords(j)[0]!=0 and intCoords(j)[0]!=(N_t-1))
			{
			for (unsigned int k=0; k<2*dim; k++) //over neighbours in both directions
				{
				double sign = pow(-1.0,k);
				int direc = (int)k/2;
				if(direc==0)
					{
               		minusDS(2*j) += real(linkMeasure/pow(dt(j),2.0))*pr(neigh(j,0,sign))-imag(linkMeasure/pow(dt(j),2.0))*pi(neigh(j,k,sign));
               		minusDS(2*j+1) += imag(linkMeasure/pow(dt(j),2.0))*pr(neigh(j,0,sign))+real(linkMeasure/pow(dt(j),2.0))*pi(neigh(j,k,sign));
					DDS.insert(2*j,2*neigh(j,0,sign)) = -real(linkMeasure/pow(dt(j),2.0));
					DDS.insert(2*j+1,2*neigh(j,0,sign))+1 = -imag(linkMeasure/pow(dt(j),2.0));
               		}
				else
					{
					minusDS(2*j) += -real(siteMeasure/pow(a,2.0))*pr(neigh(j,direc,sign))+imag(siteMeasure/pow(a,2.0))*pi(neigh(j,k,sign));
               		minusDS(2*j+1) += -imag(siteMeasure/pow(a,2.0))*pr(neigh(j,direc,sign))-real(siteMeasure/pow(a,2.0))*pi(neigh(j,k,sign));
					DDS.insert(2*j,2*neigh(j,direc,sign)) = -real(siteMeasure/pow(a,2.0));
					DDS.insert(2*j+1,2*neigh(j,direc,sign))+1 = -imag(siteMeasure/pow(a,2.0));
					}
				}
			
			minusDS(2*j) +=
			- real(siteMeasure)*( 2.0*(-(dimd-1.0)*pr(j))/pow(a,2.0) - lambda*pr(j)*(pow(pr(j),2.0)-pow(pi(j),2.0)-pow(v,2.0))/2.0 + lambda*pi(j)*(pr(j)*pi(j)-pow(v,2.0)) - epsilon/v/2.0 )
			+imag(siteMeasure)*( 2.0*(-(dimd-1.0)*pi(j))/pow(a,2.0) - lambda*pi(j)*(pow(pr(j),2.0)-pow(pi(j),2.0)-pow(v,2.0))/2.0 - lambda*pr(j)*(pr(j)*pi(j)-pow(v,2.0)) - epsilon/v/2.0 )
			- real(linkMeasure)*2.0*pr(j)/pow(a_t,2.0))
			+ imag(linkMeasure)*2.0*pi(j)/pow(a_t,2.0));
			
			minusDS(2*j+1) +=
			- imag(siteMeasure)*( 2.0*(-(dimd-1.0)*pr(j))/pow(a,2.0) - lambda*pr(j)*(pow(pr(j),2.0)-pow(pi(j),2.0)-pow(v,2.0))/2.0 + lambda*pi(j)*(pr(j)*pi(j)-pow(v,2.0)) - epsilon/v/2.0 )
			- real(siteMeasure)*( 2.0*(-(dimd-1.0)*pi(j))/pow(a,2.0) - lambda*pi(j)*(pow(pr(j),2.0)-pow(pi(j),2.0)-pow(v,2.0))/2.0 - lambda*pr(j)*(pr(j)*pi(j)-pow(v,2.0)) - epsilon/v/2.0 )
			- imag(linkMeasure)*2.0*pr(j)/pow(a_t,2.0))
			- real(linkMeasure)*2.0*pi(j)/pow(a_t,2.0));
			
			DDS.insert(2*j,2*j) = real(siteMeasure)*( 2.0*(-(dimd-1.0)/pow(a,2.0)) - lambda*(3.0*(pow(pr(j),2.0)-pow(pi(j),2.0))-pow(v,2.0))/2.0 ) -imag(siteMeasure)*lambda*3.0*pi(j)pr(j)  + real(linkMeasure)*2.0/pow(a_t,2.0);
			}
		if (intCoords(j)[0]==0)
			{
			DDS.insert(j,j) = 1.0;
			}
		if (intCoords(j)[0]==(N-1))
			{
			DDS.insert(j,j) = {0,1.0};
			DDS.insert(j,j-1) = -1/b;
			}
		} //try replacing periodic boundary conditions with phi(boundary)=v, for all boundaries.
	kinetic *= 2.0; pot_l *= 2.0; pot_e *= 2.0;//as we only calculated half a bubble in the time direction
	action = kinetic + pot_l + pot_e;

//print action early if desired
	earlyAction(print_choice, runs_count, print_run, kinetic, pot_l, pot_e, action);
//print phi, minusDS or DDS early if required
	if (print_choice != 'n' and print_choice != 'a' and runs_count==print_run)
		{
		if (print_choice=='m')
			{
			cout << left;
			//vec test;
			//test = minusDS-DDS*delta;
			for (unsigned long int j=0; j<Totdim; j++)
				{
				for (unsigned long int k=0; k<Totdim; k++)
					{
					cout << setw(7)  << setprecision(0) << fixed << DDS.coeffRef(j,k);
					}
				cout << endl;
				}
			}
		if (print_choice=='p' or print_choice == 'v')
			{
			ofstream earlyphi;
			earlyphi.open("./data/earlyphi.dat", ios::out);
			for (unsigned long int j=0; j<Totdim; j++)
				{
				for (int k=0; k<dim; k++)
					{
					earlyphi << setw(15) << coords(intCoords(j)[k]);
					}
				if (print_choice=='p')
					{
					earlyphi << setw(15) << p(j) << endl;
					}
				if (print_choice == 'v')
					{
					earlyphi << setw(15) << minusDS(j) << endl;
					}						
				}
			earlyphi.close();
			cout << "phi printed on run " << print_run << endl;
			}
		}


//solving for temp using the Newton-Raphson method
	DDS.makeCompressed();
	Eigen::SparseLU<spmat> solver;
	solver.analyzePattern(DDS);
	solver.factorize(DDS);
	delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side


//assigning values to phi
	p += delta;

//testing for where phi goes to zero - not at all correct
	for (long unsigned int j=0; j<Totdim; j++)
		{
		if( absolute(real(p(j)))<absolute(real_small) )
			{
			real_small = real(p(j));
			double real_small_loc_sqrd = 0.0;
			real_small_loc_sqrd += pow((intCoords(j)[0]-(N-1.0))/beta,2.0); //factor of 1/beta due to shorter time steps
			for (int k=1; k<dim; k++)
				{
				real_small_loc_sqrd += pow(intCoords(j)[k]-N/2.0,2.0);
				}
			small_loc = a*pow(real_small_loc_sqrd,0.5);
			}
		if( absolute(imag(p(j)))<absolute(imag_small) )
			{
			imag_small = imag(p(j));
			double imag_small_loc_sqrd = 0.0;
			imag_small_loc_sqrd += pow((intCoords(j)[0]-(N-1.0))/beta,2.0); //factor of 1/beta due to shorter time steps
			for (int k=1; k<dim; k++)
				{
				imag_small_loc_sqrd += pow(intCoords(j)[k]-N/2.0,2.0);
				}
			imag_small_loc = a*pow(imag_small_loc_sqrd,0.5);
			}
		}

bool stopper = false;
	if( runs_test<=closeness and runs_count>=min_runs) //testing on last run of newton raphson loop
		{
		log_det_DDS = solver.logAbsDeterminant(); //calculating log of determinant of DDS
		if (stopper)
			{
		for (unsigned int j=0; j<Totdim; j++) //calculating a zero mode
			{
			if (intCoords(j)[1]!=(N-1)) //doesn't currently work
				{
				dp_dx(j) = p(j+N_t)-p(j);
				}
			if (intCoords(j)[1]==(N-1))
				{
				dp_dx(j) = p(j-(N-1)*N_t)-p(j);
				}
			}
		double norm;
		norm = dp_dx.dot(dp_dx);
		vec temp_v;
		temp_v = DDS*dp_dx;
		eigen_zero = temp_v.dot(dp_dx);
		eigen_zero /= norm;
//		cout << "x-translation zero mode eigenvalue = " << eigen_zero << endl;
			}
		}

	char stop_wait;
	char print_wait;
	clock_t stop_time = clock();
	vector <int> sigma; //carries signs for whether the bubble was too large or small

	bool bool_wait = true; //set to false if you want program to stop if the looping is slow and ask the user whether or not they want to print
	convergenceQuestions (runs_count, &runs_test, min_runs, action, stop_time, &wait, &stop_wait, &print_wait, &print_choice, &print_run, kinetic, pot_l, pot_e, twaction, bool_wait, &convergence, small_loc, &gamma, sigma);

} //closing "runs" while loop

//cout << setw(8) << small_loc << setw(8) << gamma << setw(8) << action << endl;
//print_choice = 'p';
//print_run = runs_count + 1;

if (not(convergence))
	{
	runs_test = 1.0;
	}

if (convergence)
{

//stopping clock
time = clock() - time;
double realtime = time/1000000.0;

//printing to terminal
if (loop==0)
	{
	vector<string> pr_labels = {"dim","N","runs","X","action","twaction", "p=0", "time", "log|det(DDS)|", "0-mode"};
	cout_string(pr_labels);
	}
cout << left << setw(12) << dim  << setw(12) << N << setw(12) << runs_count << setw(12) << X << setw(12) << action << setw(12) << twaction << setw(12) << small_loc << setw(12) << realtime << setw(12) << log_det_DDS << setw(12) << eigen_zero << endl;

//printing action value
ofstream actionfile;
actionfile.open ("./data/action_dim.dat", ios::out | ios::app );
actionfile << left << setw(12) << dim  << setw(12) << X << setw(12) << N << setw(12) << runs_count << setw(12) << action << setw(12) << twaction << setw(12) << small_loc << setw(12) << realtime << setw(12) << log_det_DDS << setw(12) << eigen_zero << endl;
actionfile.close();

//printing output phi
string opathway = ("./data/");
string ofilename = "phi";
string oextension = (".dat");
string outfile = opathway+ofilename+to_string(loop)+oextension;
ofstream f;
f.open((outfile).c_str());
for (unsigned long int j=0; j<Totdim; j++)
	{
	intvec int_x = intCoords(j);
	vector <double> x(dim);
	f << left;
	x[0] = -box/beta + a*int_x[0]/beta; //time steps half size
	f << setw(15) << x[0];
	for (int k=1; k<dim; k++)
		{
		x[k] =  - box/2.0 + a*int_x[k];
		f << setw(15) << x[k];
		}
	f << setw(15) << p(j) << endl;	
	}
f.close();

} //closing if convergence loop
} while (not(convergence)); //closing do{}
} //closing parameter loop

return 0;
}
