#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_poly.h>
#include "parameters.h"

using namespace std;

int main()
{
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
		changeParameters (N, "N", input_p);
		}
	else if (loop_response=='y')
		{
		min_runs = 10;
		double loop_parameter = min_value + (max_value - min_value)*loop/(total_loops-1.0);
		changeParameters (loop_parameter,parameter_choice, input_p);
		}
	else if (input_p == 'b')
		{
		changeParameters (0,"null", input_p);
		}

	//starting the clock
	time = clock();
	wait = clock();

	//defining some important scalar quantities
	long double log_det_DDS;
	comp action = 2.0;
	double S_1 = 2.0*pow(mass,3)/3.0/lambda;
	double twaction = -solidAngle(dim)*epsilon*pow(R,dimd)/dimd + solidAngle(dim)*pow(R,dimd-1.0)*S_1;
	double eigen_zero;
	int alpha = 20; //gives span over which tanh is used

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = 1.0;
	int runs_count = 0;
	double runs_test = 1.0;

	//initializing phi (=p) according to choice in input_p
	vec p_e(2*Eucdim);
	vec perturb_real(Eucdim); //not sure if this is the best place to initialize the perturbations
	vec perturb_imag(Eucdim);
	vec p_neg(2*Eucdim); //negative mode of sphaleron
	
	//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
	double b_parameter = -pow(v,2.0);
	double c_parameter = epsilon/v/lambda;
	vector<double> root(3);
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());
	
	if (perturb_response=='r' or perturb_response=='i' or perturb_response=='b')
		{
		perturb_real = Eigen::VectorXd::Random(Eucdim);
		perturb_imag = Eigen::VectorXd::Random(Eucdim);
		perturb_real *= v*pow(10,-4);
		perturb_imag *= v*pow(10,-4);
		for (unsigned int j=0; j<N; j++)
			{ 
			if (input_p=='t' or input_p=='f') //explicitly 2d but generalizable
				{
				perturb_real(j*N_t) = 0;
				perturb_real((j+1)*N_t-1) = 0;
				perturb_imag(j*N_t) = 0;
				perturb_imag((j+1)*N_t-1) = 0;
				}
			else if (input_p=='b')
				{
				perturb_real((j+1)*N_t-1) += -perturb_real((j+1)*N_t-1) + perturb_real((j+1)*N_t-2);
				perturb_real(j*N_t) = 0;
				perturb_imag((j+1)*N_t-1) += -perturb_imag((j+1)*N_t-1) + perturb_imag((j+1)*N_t-2);
				perturb_imag(j*N_t) = 0;
				}
			else if (input_p == 'p')
				{
				perturb_real((j+1)*N_t-1) += - perturb_real((j+1)*N_t-1) + perturb_real((j+1)*N_t-2);
				perturb_real(j*N_t) += - perturb_real(j*N_t) + perturb_real(j*N_t+1);
				perturb_imag((j+1)*N_t-1) += -perturb_imag((j+1)*N_t-1) + perturb_imag((j+1)*N_t-2);
				perturb_imag(j*N_t) += -perturb_imag(j*N_t) + perturb_imag(j*N_t+1);
				}
			}
		}

//assigning thin wall periodic instanton solution as initial guess for phi
	for (unsigned long int j = 0; j<Eucdim; j++)
		{
		comp rho_sqrd = -pow(coord(j,0),2.0);
		comp rho1_sqrd = -pow(coord(j,0)+i*L_t/2.0,2.0);
		comp rho2_sqrd = -pow(coord(j,0)+i*L_t/2.0,2.0);
		for (int k=1; k<dim; k++)
			{
			rho_sqrd += pow(coord(j,k),2.0);
			rho1_sqrd += pow(coord(j,k)+R*cos(theta),2.0);
			rho2_sqrd += pow(coord(j,k)-R*cos(theta),2.0);
			}
		double rho = re(pow(rho_sqrd,0.5)); //should be real
		double rho1 = re(pow(rho1_sqrd,0.5));
		double rho2 = re(pow(rho2_sqrd,0.5));
		if (R<alpha/mass)
			{
			cout << "X = R*mass is too small. It should not be less than " << alpha << endl;
			}
		else
			{
			p_neg(2*j+1) = 0.0;
			p_e(2*j+1) = 0.0;
			if (input_p=='t')
				{
				p_e(2*j) = root[0];
				}
			else if (input_p=='f')
				{
				p_e(2*j) = root[2];
				}
			else if (input_p=='b')
				{
				if (rho<(R-alpha/mass))
					{
					p_e(2*j) = root[0];
					}
				else if (rho>(R+alpha/mass))
					{
					p_e(2*j) = root[2];
					}
				else
					{
					p_e(2*j) = v*tanh(mass*(rho-R)/2.0);
					p_neg(2*j) = v*pow(cosh(mass*(rho-R)/2.0),-2);
					}
				}
			else if (rho1<(R-alpha/mass) and rho2<(R-alpha/mass))
				{
				p_e(2*j) = root[0];
				}
			else if (rho1>(R+alpha/mass) or rho2>(R+alpha/mass))
				{
				p_e(2*j) = root[2];
				}
			else if (real(coord(j,1))>0) //explicitly in 2d - should be real
				{
				p_e(2*j) = v*tanh(mass*(rho1-R)/2.0);
				}
			else if (real(coord(j,1))<0) //explicitly in 2d
				{
				p_e(2*j) = v*tanh(mass*(rho2-R)/2.0);
				}
			else
				{
				p_e(2*j) = root[2]; //edge points
				}
			if (perturb_response=='r' or perturb_response=='b')
				{
				p_e(2*j) += perturb_real(j);
				}
			if (perturb_response=='i' or perturb_response=='b')
				{
				p_e(2*j+1) += perturb_imag(j);
				}
			}
		}
	
	if (input_p=='p') //fixing input phi to have zero time derivative at time boundaries
		{
		for (unsigned int j=0; j<N;j++)
			{
			p_e(2*(j*N_t+1)) += (p_e(2*j*N_t) - p_e(2*(j*N_t+1)))/2.0;
			p_e(2*j*N_t) += (-p_e(2*j*N_t) + p_e(2*(j*N_t+1)))/2.0;
			p_e(2*(j*N_t+1)+1) += (p_e(2*j*N_t+1) - p_e(2*(j*N_t+1)+1))/2.0;
			p_e(2*j*N_t+1) += (-p_e(2*j*N_t+1) + p_e(2*(j*N_t+1)+1))/2.0;
			p_e(2*((j+1)*N_t-2)) += (p_e(2*((j+1)*N_t-1)) - p_e(2*((j+1)*N_t-2)))/2.0;
			p_e(2*((j+1)*N_t-1)) += (-p_e(2*((j+1)*N_t-1)) + p_e(2*((j+1)*N_t-2)))/2.0;
			p_e(2*((j+1)*N_t-2)+1) += (p_e(2*((j+1)*N_t-1)+1) - p_e(2*((j+1)*N_t-2)+1))/2.0;
			p_e(2*((j+1)*N_t-1)+1) += (-p_e(2*((j+1)*N_t-1)+1) + p_e(2*((j+1)*N_t-2)+1))/2.0;
			}
		}
		
	//fixing norm of p_neg - may want to add a bit of p_neg to the sphaleron ('b') to make the other periodic instanton
	if (input_p == 'b')
		{
		double norm;
		norm = p_neg.dot(p_neg);
		p_neg /= norm;
		}

	//beginning newton-raphson loop, continuing loop until ratio of change of action to action is smaller than predefined "closeness".
	while (runs_test>closeness or runs_count<min_runs)
		{
			//quantities used to stop newton-raphson loop
			runs_test = abs(action - action_last)/abs(action_last);
			runs_count ++;
			action_last = action;

			// allocating memory for DS, DDS and delta
			vec minusDS(2*Eucdim);
			spMat DDS(2*Eucdim,2*Eucdim);

			int DDS_to_reserve = 2*dim+1;//number of non-zero elements per column
			DDS.reserve(Eigen::VectorXi::Constant(2*Eucdim,DDS_to_reserve));
			vec delta(2*Eucdim);  //delta=newphi-phi in Lowdim
			vec d_p_e(2*Eucdim);

			//setting action to zero, and defining quantities which split up action into its parts
			action = 0.0;
			comp kinetic = 0.0;
			comp pot_l = 0.0;
			comp pot_e = 0.0;

			//assigning values to minusDS and DDS and evaluating action
			for (unsigned long int j = 0; j < Eucdim; j++)
				{
				unsigned int time_coord = intCoords(j,N_t)[0];
				comp siteMeasure = pow(a,dim-1)*Dt(j);//for sites in time
				comp linkMeasure = pow(a,dim-1)*dt(j);//for links in time
			
				minusDS(2*j) = 0.0;//initializing to zero
				minusDS(2*j+1) = 0.0;
		
				pot_l += -siteMeasure*lambda*pow(pow(p_e(2*j)+i*p_e(2*j+1),2)-pow(v,2),2)/8.0;
				pot_e += -siteMeasure*epsilon*(p_e(2*j)+i*p_e(2*j+1)-v)/v/2.0;
				for (unsigned int k=1; k<dim; k++)
					{
					if (intCoords(j,N_t)[k]!=(N-1))
						{
						kinetic += -siteMeasure*pow(p_e(2*neigh(j,k,1,N_t))+i*p_e(2*neigh(j,k,1,N_t)+1)-p_e(2*j)-i*p_e(2*j+1),2.0)/pow(a,2)/2.0;
						}
					}
				if (time_coord==(N_t-1))
					{
					if (input_p =='p' or input_p =='b')
						{
						DDS.insert(2*j,2*j) = 1.0/b;
						DDS.insert(2*j+1,2*j+1) = 1.0;
						DDS.insert(2*j,2*(j-1)) = -1.0/b;
						}
					else if (input_p =='t' or input_p =='f')
						{
						DDS.insert(2*j,2*j) = 1.0;
						DDS.insert(2*j+1,2*j+1) = 1.0;
						}
					}
				else
					{
					kinetic += linkMeasure*pow(p_e(2*(j+1))+i*p_e(2*(j+1)+1)-p_e(2*j)-i*p_e(2*j+1),2)/pow(dt(j),2)/2.0;
					if (time_coord==0)
						{
						if (input_p =='b' or input_p =='t' or input_p =='f') //bubble boundary conditions
							{
							DDS.insert(2*j,2*j) = 1.0;
							DDS.insert(2*j+1,2*j+1) = 1.0;
							}
						else if (input_p =='p') //periodic instanton boundary conditions
							{
							DDS.insert(2*j,2*j) = -1.0/b;
							DDS.insert(2*j+1,2*j+1) = 1.0;
							DDS.insert(2*j,2*(j+1)) = 1.0/b;
							}
						}
					else
						{
						for (unsigned int k=0; k<2*dim; k++) //over neighbours in both directions
							{
							signed int sign = pow(-1,k);
							signed int deltaSign = (sign-1)/2; //zero if sign=+1 and -1 if sign=-1
							int direc = (int)k/2;
							if(direc==0)
								{
						   		minusDS(2*j) += re(pow(a,dim-1)/dt(j+deltaSign))*p_e(2*(j+sign))-im(pow(a,dim-1)/dt(j+deltaSign))*p_e(2*(j+sign)+1);
						   		minusDS(2*j+1) += re(pow(a,dim-1)/dt(j+deltaSign))*p_e(2*(j+sign)+1)+im(pow(a,dim-1)/dt(j+deltaSign))*p_e(2*(j+sign));
								DDS.insert(2*j,2*neigh(j,0,sign,N_t)) = -re(pow(a,dim-1)/dt(j+deltaSign));
								DDS.insert(2*j,2*neigh(j,0,sign,N_t)+1) = im(pow(a,dim-1)/dt(j+deltaSign));
								DDS.insert(2*j+1,2*neigh(j,0,sign,N_t)) = -im(pow(a,dim-1)/dt(j+deltaSign));
								DDS.insert(2*j+1,2*neigh(j,0,sign,N_t)+1) = -re(pow(a,dim-1)/dt(j+deltaSign));
						   		}
							else
								{
								minusDS(2*j) += -re(siteMeasure/pow(a,2))*p_e(2*neigh(j,direc,sign,N_t))+im(siteMeasure/pow(a,2))*p_e(2*neigh(j,direc,sign,N_t)+1);
								minusDS(2*j+1) += -re(siteMeasure/pow(a,2))*p_e(2*neigh(j,direc,sign,N_t)+1)-im(siteMeasure/pow(a,2))*p_e(2*neigh(j,direc,sign,N_t));
								DDS.insert(2*j,2*neigh(j,direc,sign,N_t)) = re(siteMeasure)/pow(a,2);
								DDS.insert(2*j,2*neigh(j,direc,sign,N_t)+1) = -im(siteMeasure)/pow(a,2);
								DDS.insert(2*j+1,2*neigh(j,direc,sign,N_t)) = im(siteMeasure)/pow(a,2);
								DDS.insert(2*j+1,2*neigh(j,direc,sign,N_t)+1) = re(siteMeasure)/pow(a,2);
								}				
							}
						comp temp0 = pow(a,dim-1)/dt(j) + pow(a,dim-1)/dt(j-1);
						double temp1 = -2.0*(dimd-1.0)*p_e(2*j)/pow(a,2) - lambda*p_e(2*j)*(pow(p_e(2*j),2)-3.0*pow(p_e(2*j+1),2)-pow(v,2))/2.0 - epsilon/2.0/v;
						double temp2 = -2.0*(dimd-1.0)*p_e(2*j+1)/pow(a,2) - lambda*p_e(2*j+1)*(3.0*pow(p_e(2*j),2)-pow(p_e(2*j+1),2)-pow(v,2))/2.0;
						double temp3 = -2.0*(dimd-1.0)/pow(a,2) - lambda*(3.0*pow(p_e(2*j),2)-3.0*pow(p_e(2*j+1),2)-pow(v,2))/2.0;
						double temp4 = 3.0*lambda*p_e(2*j)*p_e(2*j+1);	
					
						minusDS(2*j) += - re(siteMeasure)*temp1 + im(siteMeasure)*temp2 - re(temp0)*p_e(2*j) + im(temp0)*p_e(2*j+1);
						minusDS(2*j+1) += -im(siteMeasure)*temp1 - re(siteMeasure)*temp2 -re(temp0)*p_e(2*j+1) -im(temp0)*p_e(2*j);			
						DDS.insert(2*j,2*j) = re(siteMeasure)*temp3 + im(siteMeasure)*temp4 + re(temp0);
						DDS.insert(2*j,2*j+1) = re(siteMeasure)*temp4 - im(siteMeasure)*temp3 -im(temp0);
						DDS.insert(2*j+1,2*j) = -re(siteMeasure)*temp4 + im(siteMeasure)*temp3 + im(temp0);
						DDS.insert(2*j+1,2*j+1) = re(siteMeasure)*temp3 + im(siteMeasure)*temp4 + re(temp0);
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
					for (unsigned long int j=0; j<Eucdim; j++)
						{
						for (unsigned long int k=0; k<Eucdim; k++)
							{
							cout << setw(7)  << setprecision(0) << fixed << DDS.coeffRef(2*j,2*k);
							}
						cout << endl;
						}
					}
				else if (print_choice=='p')
					{
					earlyVector("./data/earlyphi_pic.dat",p_e);
					cout << print_choice << " printed on run " << print_run << endl;
					}
				else if (print_choice=='v')
					{
					earlyVector("./data/earlyminusDS_pic.dat",minusDS);
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
			p_e += delta;

			bool stopper = true;
			if( runs_test<=closeness and runs_count>=min_runs) //testing on last run of newton raphson loop
				{
				log_det_DDS = solver.logAbsDeterminant(); //calculating log of determinant of DDS
				if (stopper)
					{
				for (unsigned int j=0; j<Eucdim; j++)
					{
					d_p_e(2*j+1) = (p_e(2*neigh(j,1,1,N_t)+1)-p_e(2*j))/a;
					d_p_e(2*j) = (p_e(2*neigh(j,1,1,N_t))-p_e(2*j))/a;
					}
				double norm;
				norm = d_p_e.dot(d_p_e);
				vec temp_v;
				temp_v = DDS*d_p_e;
				eigen_zero = temp_v.dot(d_p_e);
				eigen_zero /= norm;
		//		cout << "x-translation zero mode eigenvalue = " << eigen_zero << endl;
					}
				}

			char stop_wait;
			char print_wait;
			clock_t stop_time = clock();

			bool bool_wait = false; //set to false if you want program to stop if the looping is slow and ask the user whether or not they want to print
			convergenceQuestions (runs_count, &runs_test, min_runs, action, stop_time, &wait, &stop_wait, &print_wait, &print_choice, &print_run, kinetic, pot_l, pot_e, twaction, bool_wait);

		} //closing "runs" while loop

	//propagating euclidean solution back in minkowskian time
	vec p_m(2*Minkdim),  vel(2*Minkdim), acc(2*Minkdim);
	
	//explicitly 2d
	for (unsigned int j=0;j<N;j++)
		{
		unsigned int k = j*N_t;
		unsigned int l = j*N_t_m;
		p_m(2*l) = p_e(2*k);
		p_m(2*l+1) = p_e(2*k+1);
		vel(2*l) = 0.0;
		vel(2*l+1) = 0.0;
		acc(2*l) = 0.0; //as the time derivatives on both sides are zero - not certain about this.
		acc(2*l+1) = 0.0;
		}
	
	for (unsigned int j=0;j<N;j++)
		{
		unsigned int l = j*N_t_m;
		p_m(2*(l+1)) = p_m(2*l) - vel(2*l)*b + acc(2*l)*pow(b,2)/2.0;
		p_m(2*(l+1)+1) = p_m(2*l+1) - vel(2*l+1)*b + acc(2*l+1)*pow(b,2)/2.0;
		}
		
	for (unsigned int j=0;j<N;j++)
		{
		unsigned int l = j*N_t_m;
		acc(2*(l+1)) += -2.0*(dimd-1.0)*p_m(2*(l+1))/pow(a,2) - lambda*p_m(2*(l+1))*(pow(p_m(2*(l+1)),2)-3.0*pow(p_m(2*(l+1)+1),2)-pow(v,2))/2.0 - epsilon/2.0/v;
		acc(2*(l+1)+1) += -2.0*(dimd-1.0)*p_m(2*(l+1)+1)/pow(a,2) - lambda*p_m(2*(l+1)+1)*(3.0*pow(p_m(2*(l+1)),2)-pow(p_m(2*(l+1)+1),2)-pow(v,2))/2.0;	
		for (unsigned int k=1; k<2*dim; k++) //over neighbours in both directions
			{
			signed int sign = pow(-1,k);
			int direc = (int)k/2;
			acc(2*(l+1)) += p_e(2*neigh(l+1,direc,sign,N_t))/pow(a,2);
			acc(2*(l+1)+1) += p_e(2*neigh(l+1,direc,sign,N_t)+1)/pow(a,2);
			}
		vel(2*(l+1)) = vel(2*l) + (acc(2*l) + acc(2*(l+1)))*b/2.0;
		vel(2*(l+1)+1) = vel(2*l+1) + (acc(2*l+1) + acc(2*(l+1)+1))*b/2.0;
		}
	
	vec p(2*Totdim); //assigning full phi vector
	for (lint j=0; j<Totdim; j++)
		{
		if (j<Minkdim)
			{
			p(2*j) = p_m(2*j);
			p(2*j+1) = p_m(2*j+1);
			}
		else
		p(2*Minkdim+2*j) = p_e(2*j);
		p(2*Minkdim+2*j+1) = p_e(2*j+1);
		}

	//stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;

	//printing to terminal
	if (loop==0)
		{
		vector<string> short_labels = {"dim","N","X","runs", "time"};
		vector<string> long_labels = {"re(action)", "im(action)", "log|det(DDS)|", "0-mode"};
		coutStringShort(short_labels);
		coutStringLong(long_labels);
		}
	cout << left << setw(8) << dim  << setw(8) << N << setw(8) << X << setw(8) << runs_count << setw(8) << realtime << setw(20) << re(action) << setw(20) << im(action) << setw(20) << log_det_DDS << setw(20) << eigen_zero << endl;

	//printing action value
	ofstream actionfile;
	actionfile.open ("./data/action_pic.dat", ios::out | ios::app );
	actionfile << left << setw(8) << dim  << setw(8) << N << setw(8) << X << setw(8) << runs_count << setw(8) << realtime << setw(20) << re(action) << setw(20) << im(action) << setw(20) << log_det_DDS << setw(20) << eigen_zero << endl;
	actionfile.close();

	//printing output phi
	string opathway = ("./data/");
	string ofilename = "phi_pic";
	string oextension = (".dat");
	string outfile = opathway+ofilename+to_string(loop)+oextension;
	ofstream f;
	f.open((outfile).c_str());
	for (unsigned long int j=0; j<Eucdim; j++)
		{
		f << left;
		for (int r=0; r<dim; r++)
			{
			f << setw(15) << re(coord(j,r)) << setw(15) << im(coord(j,r));
			}
		f << setw(15) << p(2*j) << setw(15) << p(2*j+1) << endl;	
		}
	f.close();

} //closing parameter loop

return 0;
}