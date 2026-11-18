// Calculates mean squared displacement of a trajectory from the Brownian Motion simulation

// Libraries
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

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
    // Read in trajectory from .xyz file produced by LAMMPS
    std::vector<double> x_trajectory;
    x_trajectory.reserve(100000001);
    std::ifstream file("./unwrapped_trajectory_3d.xyz");
    std::string line;
    int counter = 1;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double first_num;
        
        if (iss >> first_num) {
            if (counter % 10 == 0) {
                x_trajectory.push_back(first_num);
            }
        }
        counter=counter+1;
    }

    // Calculate mean squared displacement
    int max_lag {10000000};
    int num_timesteps {static_cast<int>(size(x_trajectory))};
    std::vector<double> msd {MeanSquaredDisplacement(x_trajectory,max_lag,num_timesteps)};

    // Write mean squared displacements to a file
    std::ofstream msd_file("./msd_he_ar_3d.bin", std::ios::binary);
    msd_file.write(reinterpret_cast<const char*>(msd.data()),msd.size()*sizeof(double));

    return 0;
}
