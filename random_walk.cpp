// Calculates random walk and mean squared displacement from that random walk

// Libraries
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <random>

// Global constants
const double boltzmann {1.380649e-23};

// Function for Brownian Random Walk
std::vector<double> RandomWalk(double x0, double x1, double deltat, double D, double T, double mass, int num_timesteps, 
    std::mt19937 generator, std::normal_distribution<double> d){
    
    // Set array to store trajectory
    std::vector<double> x_trajectory (num_timesteps);
    x_trajectory[0] = x0;
    x_trajectory[1] = x1;

    // Recursively calculate positions and assign to appropriate element of x_trajectory
    for(int i {2}; i<num_timesteps; i=i+1){
        double noise {d(generator)};
        double ratio {deltat*boltzmann*T/(D*mass)};
        x_trajectory[i] = x_trajectory[i-1]*((2 + ratio)/(1+ratio)) - 1/(1+ratio)*x_trajectory[i-2] + boltzmann*T*sqrt(2/D)/(mass*(1+ratio))*pow(deltat,1.5)*noise;
    }
    return x_trajectory;
}

// Function to calculate mean squared displacements 
std::vector<double> MeanSquaredDisplacement(std::vector<double> x_trajectory, int max_lag, int num_timesteps){

    // Vector to record mean squared displacements
    std::vector<double> mean_squared_displacements (max_lag);

    // Loop through the lags
    for(int lag{1}; lag <= max_lag; lag=lag+100){
        double msd {0.0};

        // Loop through the timesteps to sum up the squared displacements at a particular lag
        for(int timestep{0}; timestep<num_timesteps-lag; timestep=timestep+1){
            msd = msd + (x_trajectory[timestep+lag] - x_trajectory[timestep])*(x_trajectory[timestep+lag] - x_trajectory[timestep]);
        }
        // Divide by the number of squared displacements to get the mean
        msd = msd/(num_timesteps-lag);

        mean_squared_displacements[lag-1] = msd;
    }

    return mean_squared_displacements;
}

int main(){

    // Lennard Jones constants
    double epsilon {9.8725*boltzmann};
    double sigma {252.38*1e-12};
    double mass_he {6.64e-27};
    double mass_ar {6.63e-26};

    // Set simulation parameters
    double x0 {0};
    double x1 {0};
    double D {6.4e-5};
    double deltat {0.001*sqrt(mass_he*sigma*sigma/epsilon)};
    double T {298};
    int num_timesteps {100000000};

    // Create PRNG for the random noise
    std::mt19937 generator{};
    std::normal_distribution<double> d{0.0,1.0};

    // Calculate random walk
    std::vector<double> x_trajectory {RandomWalk(x0,x1,deltat,D,T,mass_ar,num_timesteps,generator,d)};

    // Calculate mean squared displacement from the trajectory
    int max_lag {10000000};
    std::vector <double> msd {MeanSquaredDisplacement(x_trajectory,max_lag,num_timesteps)};

    // Write x trajectory to file for post-processing
    std::ofstream xfile("./x_trajectory_random_walk_2.bin", std::ios::binary);
    xfile.write(reinterpret_cast<const char*>(x_trajectory.data()),x_trajectory.size()*sizeof(double));

    // Write mean squared displacements to a file
    std::ofstream msd_file("./msd_random_walk_2b.bin", std::ios::binary);
    msd_file.write(reinterpret_cast<const char*>(msd.data()),msd.size()*sizeof(double));
}
