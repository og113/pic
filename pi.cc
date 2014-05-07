#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_poly.h>
#include "parameters1.h"

using namespace std;
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
char input_p;
char perturb_response;
char loop_response;
string parameter_choice;
int min_value;
int max_value;
int total_loops;
char print_choice;
int print_run;
int min_runs;

//asking user questions
askQuestions (&input_p, &perturb_response, &loop_response, &parameter_choice, &min_value, &max_value, &total_loops, &print_choice, &print_run, &min_runs);

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
	change_parameters (N, "N", input_p);
	}
else if (loop_response=='y')
	{
	min_runs = 10;
	double loop_parameter = min_value + (max_value - min_value)*loop/(total_loops-1.0);
	change_parameters (loop_parameter,parameter_choice, input_p);
	}
else if (input_p == 'b')
	{
	change_parameters (0,"null", input_p);
	}

//starting the clock
time = clock();
wait = clock();

//defining some important scalar quantities
long double log_det_DDS;
double action = 2.0;
double S_1 = epsilon*R/(dimd-1);
double twaction = -solid_angle*epsilon*pow(R,dimd)/dimd + solid_angle*pow(R,dimd-1.0)*S_1;
double small_loc = 0;
double eigen_zero;
//alpha gives the range over which the thin wall solution tanh is used
double alpha = 20.0;
//gamma gives an extra bit of size to the bubble
double gamma = 0.0;

//defining some quantities used to stop the Newton-Raphson loop when action stops varying
double action_last = 1.0;
int runs_count = 0;
double runs_test = 1.0;

bool convergence = false;
do {
//initializing phi (=p)
vec p_e(Eucdim);
vec perturb(Eucdim);

//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
double b_parameter = -pow(v,2.0);
double c_parameter = epsilon/v/lambda;
vector<double> root(3);
gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
sort(root.begin(),root.end());

if (perturb_response=='y')
		{
		perturb = Eigen::VectorXd::Random(Eucdim);
		perturb *= v*pow(10,-4);
		for (unsigned int j=0; j<N; j++)
			{
			if (input_p=='t' or input_p=='f') //explicitly 2d but generalizable
				{
				perturb(j*N_te) = 0;
				perturb((j+1)*N_te-1) = 0;
				}
			else if (input_p=='b')
				{
				perturb((j+1)*N_te-1) += -perturb((j+1)*N_te-1) + perturb((j+1)*N_te-2);
				perturb(j*N_te) = 0;
				}
			else if (input_p == 'p')
				{
				perturb((j+1)*N_te-1) += -perturb((j+1)*N_te-1) + perturb((j+1)*N_te-2);
				perturb(j*N_te) += -perturb(j*N_te) + perturb(j*N_te+1);
				}
			}
		}


//assigning thin wall periodic instanton solution as initial guess for phi
for (unsigned long int i = 0; i<Eucdim; i++)
	{
	double bubble_size = R * (1.0 + gamma);
	double rho_sqrd = pow(coord(i,0),2.0);
	double rho1_sqrd = pow(coord(i,0)+L_te/2.0,2.0);
	double rho2_sqrd = pow(coord(i,0)+L_te/2.0,2.0);
	for (int k=1; k<dim; k++)
		{
		rho_sqrd += pow(coord(i,k),2.0);
		rho1_sqrd += pow(coord(i,k)+R*cos(theta),2.0);
		rho2_sqrd += pow(coord(i,k)-R*cos(theta),2.0);
		}
	double rho = pow(rho_sqrd,0.5);
	double rho1 = pow(rho1_sqrd,0.5);
	double rho2 = pow(rho2_sqrd,0.5);
	if (R<alpha/mass)
		{
		cout << "X = R*mass is too small. It should not be less than " << alpha << endl;
		}
	else
		{
		if (input_p=='t')
			{
			p_e(i) = root[0];
			}
		else if (input_p=='f')
			{
			p_e(i) = root[2];
			}
		else if (input_p=='b')
			{
			if (rho<(bubble_size-alpha/mass))
				{
				p_e(i) = root[0];
				}
			else if (rho>(bubble_size+alpha/mass))
				{
				p_e(i) = root[2];
				}
			else
				{
				p_e(i) = v*tanh(mass*(rho-bubble_size)/2.0);
				}
			}
		else if (rho1<(bubble_size-alpha/mass) and rho2<(bubble_size-alpha/mass))
			{
			p_e(i) = root[0];
			}
		else if (rho1>(bubble_size+alpha/mass) or rho2>(bubble_size+alpha/mass))
			{
			p_e(i) = root[2];
			}
		else if (coord(i,1)>0) //explicitly in 2d
			{
			p_e(i) = v*tanh(mass*(rho1-bubble_size)/2.0);
			}
		else if (coord(i,1)<0) //explicitly in 2d
			{
			p_e(i) = v*tanh(mass*(rho2-bubble_size)/2.0);
			}
		else
			{
			p_e(i) = root[2]; //edge points
			}
		if (perturb_response=='y')
				{
				p_e(i) += perturb(i);
				}
		}
	}
	
if (input_p=='p') //fixing input phi to have zero time derivative at time boundaries
	{
	for (unsigned int j=0; j<N;j++)
		{
		p_e(j*N_te+1) += (p_e(j*N_te) - p_e(j*N_te+1))/2.0;
		p_e(j*N_te) += (-p_e(j*N_te) + p_e(j*N_te+1))/2.0;
		p_e((j+1)*N_te-2) += (p_e((j+1)*N_te-1) - p_e((j+1)*N_te-2))/2.0;
		p_e((j+1)*N_te-1) += (-p_e((j+1)*N_te-1) + p_e((j+1)*N_te-2))/2.0;
		}
	}
	
//beginning newton-raphson loop, continuing loop until ratio of change of action to action is smaller than predefined "closeness".
while (runs_test>closeness or runs_count<min_runs)
{
//quantities used to stop newton-raphson loop
runs_test = absolute(action - action_last)/absolute(action_last);
runs_count ++;
action_last = action;

// allocating memory for DS, DDS and delta
vec minusDS(Eucdim);
spMat DDS(Eucdim,Eucdim);
int DDS_to_reserve = 2*dim+1;//number of non-zero elements per column
DDS.reserve(Eigen::VectorXi::Constant(Eucdim,DDS_to_reserve));
vec delta(Eucdim);  //delta=newphi-phi in Lowdim
vec d_p_e(Eucdim);

//setting action to zero, and defining quantities which split up action into its parts
action = 0.0;
double kinetic = 0.0;
double pot_l = 0.0;
double pot_e = 0.0;
double small = p_e(0);
//assigning values to minusDS and DDS and evaluating action
	for (unsigned long int i = 0; i < Eucdim; i++)
		{
		unsigned int time_coord = intCoords(i)[0];
		double measure = pow(a,dimd-1.0)*b;
		minusDS(i) = 0.0;//initializing to zero
		if (time_coord==(N_te-1))
			{
			pot_l += 0.5*measure*lambda*pow(pow(p_e(i),2)-pow(v,2),2.0)/8.0;
			pot_e += 0.5*measure*epsilon*(p_e(i)-v)/v/2.0;
			if (input_p =='b' or input_p =='p')
					{
					DDS.insert(i,i) = 1.0/b;
					DDS.insert(i,i-1) = -1.0/b;
					}
				else if (input_p =='t' or input_p =='f')
					{
					DDS.insert(i,i) = 1.0;
					}
			}
		else 
			{
			pot_l += measure*lambda*pow(pow(p_e(i),2)-pow(v,2),2.0)/8.0;
			pot_e += measure*epsilon*(p_e(i)-v)/v/2.0;
	      		for (unsigned int j=0; j<dim; j++) //just over neighbours in positive directions
	           		{
					if(j==0)
						{			
						kinetic += measure*pow(p_e(neigh(i,0,1))-p_e(i),2.0)/pow(b,2.0)/2.0;
						}
					else
						{
						kinetic += measure*pow(p_e(neigh(i,j,1))-p_e(i),2.0)/pow(a,2.0)/2.0;
						}
	           		}
	        if (time_coord==0)
				{
				pot_l += -0.5*measure*lambda*pow(pow(p_e(i),2)-pow(v,2),2)/8.0;
				pot_e += -0.5*measure*epsilon*(p_e(i)-v)/v/2.0;
				if (input_p =='b' or input_p =='t' or input_p =='f')
					{
					DDS.insert(i,i) = 1.0;
					}
				else if (input_p =='p')
					{
					DDS.insert(i,i) = -1.0/b;
					DDS.insert(i,i+1) = 1.0/b;
					}
				}
			else
				{
				for (unsigned int k=0; k<2*dim; k++) //over neighbours in both directions
					{
					signed int sign = pow(-1,k);
					int direc = (int)k/2;
					if(direc==0)
						{
		           		minusDS(i) += measure*p_e(i+sign)/pow(b,2.0);
						DDS.insert(i,i+sign) = -measure/pow(b,2.0);
		           		}
					else
						{
						minusDS(i) += measure*p_e(neigh(i,direc,sign))/pow(a,2.0);
						DDS.insert(i,neigh(i,direc,sign)) = -measure/pow(a,2.0);
						}				
					}
				minusDS(i) += - measure*( (2.0*(dimd-1.0)/pow(a,2.0) + 2.0/pow(b,2.0))*p_e(i) + lambda*p_e(i)*(pow(p_e(i),2.0)-pow(v,2.0))/2.0 + epsilon/v/2.0 );
				DDS.insert(i,i) = measure*( 2.0*(dimd-1.0)/pow(a,2.0) + 2.0/pow(b,2.0) + lambda*(3.0*pow(p_e(i),2.0)-pow(v,2.0))/2.0 );
				}
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
			for (unsigned long int i=0; i<Eucdim; i++)
				{
				for (unsigned long int j=0; j<Eucdim; j++)
					{
					cout << setw(7)  << setprecision(0) << fixed << DDS.coeffRef(i,j);
					}
				cout << endl;
				}
			}
		else if (print_choice=='p' or print_choice == 'v')
			{
			ofstream earlyphi;
			earlyphi.open("./data/earlyphi.dat", ios::out);
			for (unsigned long int i=0; i<Eucdim; i++)
				{
				earlyphi << left;
				for (int r=0; r<dim; r++)
					{
					earlyphi << setw(15) << coord(i,r);
					}
				if (print_choice=='p')
					{
					earlyphi << setw(15) << p_e(i) << endl;
					}
				if (print_choice == 'v')
					{
					earlyphi << setw(15) << minusDS(i) << endl;
					}						
				}
			earlyphi.close();
			cout << print_choice << " printed on run " << print_run << endl;
			}
		}


//solving for temp using the Newton-Raphson method
	DDS.makeCompressed();
	Eigen::SparseLU<spMat> solver;
	solver.analyzePattern(DDS);
	solver.factorize(DDS);
	delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side


//assigning values to phi
	p += delta;
	
	bool smallLoc = false;

//testing for where phi goes to zero
	if (smallLoc)
		{
		for (long unsigned int i=0; i<Eucdim; i++)
			{
			if( absolute(p_e(i))<absolute(small) )
				{
				small = p_e(i);
				double small_loc_sqrd = 0.0;
				small_loc_sqrd += pow((coord(i,0)-b*(N-1.0)),2.0); //this is completely wrong
				for (int j=1; j<dim; j++)
					{
					small_loc_sqrd += pow(coord(i,j)-N/2.0,2.0);
					}
				small_loc = a*pow(small_loc_sqrd,0.5);
				}
			}
		}

bool stopper = true;
	if( runs_test<=closeness and runs_count>=min_runs) //testing on last run of newton raphson loop
		{
		log_det_DDS = solver.logAbsDeterminant(); //calculating log of determinant of DDS
		if (stopper)
			{
		for (unsigned int i=0; i<Eucdim; i++)
			{
			d_p_e(i) = p_e(neigh(i,1,1))-p_e(i);
			}
		double norm;
		norm = d_p.dot(d_p);
		vec temp_v;
		temp_v = DDS*d_p;
		eigen_zero = temp_v.dot(d_p);
		eigen_zero /= norm;
//		cout << "x-translation zero mode eigenvalue = " << eigen_zero << endl;
			}
		}

	char stop_wait;
	char print_wait;
	clock_t stop_time = clock();
	vector <int> sigma; //carries signs for whether the bubble was too large or small

	bool bool_wait = false; //set to false if you want program to stop if the looping is slow and ask the user whether or not they want to print
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
actionfile.open ("./data/action_pi.dat", ios::out | ios::app );
actionfile << left << setw(12) << dim  << setw(12) << X << setw(12) << N << setw(12) << runs_count << setw(12) << action << setw(12) << twaction << setw(12) << small_loc << setw(12) << realtime << setw(12) << log_det_DDS << setw(12) << eigen_zero << endl;
actionfile.close();

//printing output phi
string opathway = ("./data/");
string ofilename = "phi";
string oextension = (".dat");
string outfile = opathway+ofilename+to_string(loop)+oextension;
ofstream f;
f.open((outfile).c_str());
for (unsigned long int i=0; i<Eucdim; i++)
	{
	f << left;
	for (int r=0; r<dim; r++)
		{
		f << setw(15) << coord(i,r);
		}
	f << setw(15) << p_e(i) << endl;	
	}
f.close();

} //closing if convergence loop
} while (not(convergence)); //closing do{}
} //closing parameter loop

return 0;
}
