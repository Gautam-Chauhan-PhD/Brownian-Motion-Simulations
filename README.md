# Brownian-Motion-Simulations
A repository containing simulations of Brownian Motion. The system we model consists of an Argon atom diffusing through Helium gas, and undergoing Brownian Motion due to collisions with the Helium atoms. We model this system using molecular dynamics with a Lennard-Jones force to model the He-Ar interactions. We then calculate the mean squared displacement of the Argon atom in the x direction for various lags and plot this against the lag. The diffusion coefficient can be estimated from the linear portion of the graph using the equation:

$MSD = 2D \Delta t$

Finally, we compare the diffusion coefficient to the value from experiments and compare the mean squared displacement graph to that obtained with Langevin dynamics using the experimental diffusion coefficient.

The main purpose of this mini-project was to learn the following programming tools:
- The molecular dynamics simulation was implemented using LAMMPS
- The Langevin dynamics simulation was implemented from scratch using C++
- Post-processing of data and plotting was done using python in Jupyter Notebook
