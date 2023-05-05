//
//  main.cpp
//  chain
//
//  Created by Aleksei Sokolov on 03.08.18.
//  Copyright © 2018 Aleksei Sokolov. All rights reserved.
//
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
using namespace std;

//simulation parameters
int N = 201;
int Ntime = 5000;
int Nrealisations = 1000;
double C = 200;
double a = 1.0;
double m = 1.0;
double k_b = 1.0;
double dt = 0.02*M_PI*sqrt(m/C);


//open output file
int output_interval = 100;
ofstream file_data;

vector<double> ddu (N);
vector<double> du (N);
vector<double> u (N);

double pair_force(double z, double a, double C) {
    double alpha = 500.0;
    double beta = 0.0;
    double force = C*(a-z) + alpha*pow((a-z), 2.0) + beta*pow((a-z), 3.0);
    return force;
}

double initial_displacement (double n, double a) {
	double u_0 = 0;
	//u_0 = exp(-pow((n-int(N/2)), 2)/50)*sin(1*n); 
    u_0 =sin(2*3.14/(N-1)*n);
    return u_0;
}

double initial_velocities (double n, double a) {
	double u_0 = 0;
	u_0 = exp(-pow((n-int(N/2)), 2)/50)*cos(1*n); 
    return u_0;
}

int main(int argc, const char * argv[]) {

    for (int j = 0; j<N; j++) {
		u[j] = initial_displacement(j,a);
		du[j] = initial_velocities(j,a);
	}

    //loop
    for (int t=0; t<Ntime; t++) {
        if(t % output_interval == 0) {
	        string filename;		
			filename = "data/" + to_string(t) + ".dump";
			file_data.open (filename);
            for(int j=0; j<N; j++) file_data << j << " " << fixed << u[j] <<"\n";
			file_data.close();
		}
        // calculate forces and next step over realisations
		ddu[0] = - pair_force(u[1] - u[0], a, C) + pair_force(u[0] - u[N-1], a, C) ;  //periodic left
		
		//forces
		for(int j=1; j < N-1; j++) {
			ddu[j] = - pair_force(u[j+1]- u[j], a, C)  + pair_force(u[j] - u[j-1], a, C) ;
		}
		
		ddu[N-1] = - pair_force(u[0] - u[N-1], a, C)  + pair_force(u[N-1] - u[N-2], a, C); //periodic right
		
		//integration
		for(int j=0; j<N; j++){
			du[j] += dt*ddu[j];
			u[j] += dt*du[j];
        }
    }

    
    return 0;
}

