//
//

#ifndef MINIMD_SIMULATION_HPP
#define MINIMD_SIMULATION_HPP

#include "eigen3/Eigen/Dense"

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

class Simulation {
private:
    // Constants
    const unsigned runTime_;
    const unsigned nAtoms_;
    const double mass_ = 1;
    const double kb_ = 1;
    const double rc_ = 3.5294;
    const double volume_ = 5.32186*5.32186*5.32186;  // In reduced len^3

    // Simulation box geometry
    Eigen::RowVector3d xVec_ = {5.32186,  0,  0};
    Eigen::RowVector3d yVec_ = {0,  5.32186,  0};
    Eigen::RowVector3d zVec_ = {0,  0,  5.32186};

    // These quantities will change throughout the simulation
    double temp_;
    double pressure_;
    double potential_energy_ = 0;
    double kinetic_energy_ = 0;
    double density_ = 0;
    double virial_ = 0;
    int time_ = 0;


public:
    Eigen::Matrix<double, Eigen::Dynamic, 3> vel_;
    Eigen::Matrix<double, Eigen::Dynamic, 3> pos_;


    // Getters
    int nAtoms() {return nAtoms_;}
    double dim() {return xVec_(0);}
    double rc() {return rc_;}
    double mass() {return mass_;}
    double dens() {return density_;}

    // Setters
    void setKE(double KE) {kinetic_energy_ = KE;}
    void setPE(double PE) {potential_energy_ = PE;}
    void setVirial(double vir) {virial_ = vir;}
    void setPress(double p) {pressure_ = p;}
    void setTemp(double temp) {temp_ = temp;}
    void incTime(){time_++;}

    // Constructor
    Simulation(unsigned nAtoms, double temp, unsigned runTime);

    // Outputs
    void print_pos();
    void print_stat();
};


// Integrator for the simulation
class VelocityVerletIntegrator {
private:
    Eigen::Matrix<double, Eigen::Dynamic, 3> forces;
    Simulation& sim_;
    int nAtoms = 0;
    double dim = 0;
    double rc = 0;
    double E_pot = 0;
    double E_kin = 0;
    double virial = 0;
    double press = 0;
    double temp = 0;

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

    void VelocityVerletStep(double ts);

    void SetValues();

};

#endif //MINIMD_SIMULATION_HPP
