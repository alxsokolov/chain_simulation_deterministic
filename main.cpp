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
int N = 501;
int Ntime = 50000;
double C = 200;
double a = 1.0;
double m = 1.0;
double k_b = 1.0;
double sigma = 50;
double dt = 0.02*M_PI*sqrt(m/C);
double sound_vel = a*sqrt(C/m);
double k0 = 2*M_PI/(N-1)*40;
double k1 = 2*M_PI/(N-1)*20;
double wp_dist = 50;


//open output file
int output_interval = 5000;
ofstream file_data;

vector<double> ddu (N);
vector<double> du (N);
vector<double> u (N);

double pair_force(double z, double a, double C) {
    double alpha = 1.0;
    double beta = 0.0;
    double force = C*(a-z) + alpha*pow((a-z), 2.0) + beta*pow((a-z), 3.0);
    return force;
}

double initial_displacement (double n, double a) {
	double u_0 = 0;
	double u_1 = 0;
	u_0 = exp(-pow((n-(N-1)/2+wp_dist), 2)/2/pow(sigma,2))*sin(k0*n); 
	u_1 = exp(-pow((n-(N-1)/2-wp_dist), 2)/2/pow(sigma,2))*sin(k1*n); 
    return u_1;
}

double initial_velocities (double n, double a) {
	double u_0 = 0;
	double u_1 = 0;
	u_0 = exp(-pow((n-(N-1)/2+wp_dist), 2)/2/pow(sigma,2))*k0*cos(k0*n) - pow(sigma, -2)*exp(-pow((n-(N-1)/2)+wp_dist, 2)/2/pow(sigma,2))*(n-(N-1)/2+wp_dist)*sin(k0*n); 
	u_1 = exp(-pow((n-(N-1)/2)-wp_dist, 2)/2/pow(sigma,2))*k1*cos(k1*n) - pow(sigma, -2)*exp(-pow((n-(N-1)/2) -wp_dist, 2)/2/pow(sigma,2))*(n-(N-1)/2-wp_dist)*sin(k1*n); 
	u_0 *= -sound_vel;
	u_1 *= sound_vel;
    return u_1;
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

