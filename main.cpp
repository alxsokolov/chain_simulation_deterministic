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
int N = 200;
int Ntime = 1000;
int Nrealisations = 1000;
double C = 200;
double a = 1.0;
double m = 1.0;
double k_b = 1.0;
double dt = 0.02*M_PI*sqrt(m/C);
double bath_width = 5;

//Langeven thermostat parameters
double Thot = 1;
double Tcold = 0.5;
double lambda_hot = 10;
double lambda_cold = 1;

//open output file
int output_interval = (int) 10/(dt*sqrt(C/m));
ofstream file_data;

vector<vector<double>> ddu (Nrealisations, vector<double> (N));
vector<vector<double>> du (Nrealisations, vector<double> (N));
vector<vector<double>> u (Nrealisations, vector<double> (N));
vector<double> T (N);
vector<double> h (N);

double pair_force(double z, double a, double C) {
    double alpha = 0.0;
    double beta = 0.0;
    double force = C*(a-z) + alpha*pow((a-z), 2.0) + beta*pow((a-z), 3.0);
    return force;
}

double initial_temp_perturbation (double x, int N, double a) {
    double T_0 = 0;
    if(x<a*N/2) T_0 = 1;
    return T_0;
}

double initial_temp_perturbation_uniform (double x, int N, double a, double T) {
    return T;
}

int main(int argc, const char * argv[]) {
    
    file_data.open ("data.txt");
    
    //initial conditions stochastic
    for (int i=0; i < Nrealisations; i++){
        double width = 100;
        double mean = 0;
        double sigma = 1;
        random_device rd;
        default_random_engine generator(rd()); //seed
        normal_distribution<double> nd(mean, sigma);
        //for (int j = N/2 - width/2; j < N/2 + width/2; j++) du[i][j] = nd(generator);
        //thermal noise over the whole chain
        //for (int j = 1; j < N-1; j++) du[i][j] = 10*nd(generator);
        for (int j = 0; j<N; j++) du[i][j] = sqrt(2*k_b*initial_temp_perturbation_uniform(a*j,N,a,0)/m) * nd(generator);
    }

    //loop
    for (int t=0; t<Ntime; t++) {
        
        //output to file
       // if((int) t*dt*sqrt(C/m) >5 && (int) t*dt*sqrt(C/m) <5.05)
       //     for(int j=0; j<N; j++) file_data << t*dt*sqrt(C/m) << " " << j << " " << fixed << T[j] << " " << h[j] <<"\n";
        
        if(t % output_interval == 0)
            for(int j=0; j<N; j++) file_data << t*dt*sqrt(C/m) << " " << j << " " << fixed << T[j] << " " << h[j] <<"\n";
        // averaged quntities to zero
        for(int j=0; j<N; j++) {T[j] = 0; h[j] = 0;}
        // calculate forces and next step over realisations
        for (int i = 0; i<Nrealisations; i++) {
            //boundary conditions left
            
            //ddu[i][0] =  - pair_force(u[i][1] - u[i][0]); //free left
            ddu[i][0] = - pair_force(u[i][1] - u[i][0], a, C) + pair_force(u[i][0] - u[i][N-1], a, C) ;  //periodic left
            //ddu[i][0]= 0; //fixed left
            
            //forces
            for(int j=1; j < N-1; j++) {
                ddu[i][j] = - pair_force(u[i][j+1]- u[i][j], a, C)  + pair_force(u[i][j] - u[i][j-1], a, C) ;
            }
            
            //boundary conditions right
            //ddu[i][N-1] = pair_force(u[i][N-1] - u[i][N-2]); //free right
            ddu[i][N-1] = - pair_force(u[i][0] - u[i][N-1], a, C)  + pair_force(u[i][N-1] - u[i][N-2], a, C); //periodic right
            //ddu[i][N-1] = 0; //fixed right
            
            //add Langeven thermostats
            double n_left = (int) N/4;
            double n_right = (int) 3*N/4;
            for(int h = -(int)bath_width/2; h < (int)bath_width/2; h++)
            {
                random_device lngvn_hot;
                default_random_engine gen_lngvn_hot(lngvn_hot()); //seed
                random_device lngvn_cold;
                default_random_engine gen_lngvn_cold(lngvn_cold()); //seed
                normal_distribution<double> Xi_hot(0, 2*lambda_hot*k_b*Thot);
                normal_distribution<double> Xi_cold(0, 2*lambda_cold*k_b*Tcold);
                ddu[i][n_left+h] += Xi_hot(gen_lngvn_hot) - lambda_hot*du[i][n_left+h];
                //ddu[i][n_right+h] += Xi_cold(gen_lngvn_cold) - lambda_cold*du[i][n_right+h];
            }
            
            //integration
            for(int j=0; j<N; j++){
                du[i][j] += dt*ddu[i][j];
                u[i][j] += dt*du[i][j];
            }
        }
        //averaging
        for(int j=0; j<N; j++)  for(int i=0; i<Nrealisations; i++) T[j] += m / k_b * pow(du[i][j],2);
        for(int j=0; j<N-1; j++)  for(int i=0; i<Nrealisations; i++) h[j] += 0.5*a*(du[i][j+1] + du[i][j])*pair_force(u[i][j+1]- u[i][j], a, C);
        for(int j=0; j<N; j++) {
            T[j] = T[j]/Nrealisations;
            h[j] = h[j]/Nrealisations;
        }
        //percents done
        if(t % (Ntime / 10) == 0) cout<< (double) t/(double) Ntime * 100 << "% \n";
    }

    file_data.close();
    
    
    std::cout << "Hello, World!\n";
    return 0;
}

