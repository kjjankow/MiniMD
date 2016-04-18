//
//

#ifndef MINIMD_SIMULATION_HPP
#define MINIMD_SIMULATION_HPP

#include "eigen3/Eigen/Dense"
#include "helper.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

class Simulation {
private:
    // Constants
    const unsigned runTime_;
    const unsigned nAtoms_;
    const double mass_ = 39.948;            // Hardcode Argon's atomic mass
    const double kb_ = 1.380648e-23;        // Boltzmann's const in J/K
    const double rc_ = 12.0;                // Cutoff radius in Ang
    const double volume_ = 30.0*30.0*30.0;  // In Ang^3. Traditionally we would take the scalar triple product
                                            // of the three dimension vectors, but we assume a cubic box here.

    // Simulation box geometry
    Eigen::RowVector3d xVec_ = {30.0,    0,    0};
    Eigen::RowVector3d yVec_ = {0   , 30.0,    0};
    Eigen::RowVector3d zVec_ = {0   ,    0, 30.0};

    // These quantities will change throughout the simulation
    double temp_;
    double pressure_;


public:
    Eigen::Matrix<double, Eigen::Dynamic, 3> vel_;
    Eigen::Matrix<double, Eigen::Dynamic, 3> pos_;
    // Constructor
    Simulation(unsigned nAtoms, double temp, unsigned runTime);

    // Getters
    int nAtoms() {return nAtoms_;}
    double dim() {return xVec_(1);}
    double rc() {return rc_;}

    // Setters

    // Outputs
    void print_pos();
};


// Integrator for the simulation
class VelocityVerletIntegrator {
private:
    Eigen::Matrix<double, Eigen::Dynamic, 3> forces;
    Eigen::Matrix<double, Eigen::Dynamic, 3> accel;
    Simulation& sim_;
    int nAtoms = 0;
    double E_pot = 0;
    double virial = 0;

public:
    // Constructor
    VelocityVerletIntegrator(Simulation& sim);

    // Propagation functions

    /*
     * Calculates the forces due to Lenanrd-Jones (LJ) interactions between
     * all atom pairs from the simulation pointed at in the integrator. Sigma is an experimentally
     * determined scaling parameter which describes the equilibrium separation distance between two
     * atoms. Eps is an experimentally determined scaling parameter which gives the energy of this
     * interaction at at distance of sigma. Uses the standard minimum image convention.
     */
    void Forces();

    /*
     * Wraps the position of any atom that leaves the simulation box back into the opposite side
     * to enforce periodic boundary conditions.
     */
    void WrapPos();

};

#endif //MINIMD_SIMULATION_HPP
