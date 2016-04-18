//
// Created by Tom Parsons on 4/14/16.
//

#ifndef MINIMD_HELPER_HPP
#define MINIMD_HELPER_HPP

#include <cmath>
#include "eigen3/Eigen/Dense"

/*
 * Calculates the distance vector between two atoms under periodic boundary conditions
 * according to the minimum image criterion. Assumes a cubic box.
 */
Eigen::Vector3d PeriodicDistanceVec(Eigen::Vector3d x_pos, Eigen::Vector3d y_pos, double dim);


/*
 * Calculates the force due to a Lenanrd-Jones (LJ) interaction between
 * a pair of atoms from the distance vector between them. Sigma is an experimentally
 * determined scaling parameter which describes the equilibrium separation distance between two
 * atoms. Eps is an experimentally determined scaling parameter which gives the energy of this
 * interaction at at distance of sigma.
 *
 * \param dist_vec Distance between an atom pair in Angstroms
 * \param rc Cutoff radius in Angstroms
 * \param sigma LJ parameter in Angstroms
 * \param eps LJ parameter in
 */
Eigen::Vector3d LJForce(Eigen::Vector3d dist_vec, double rc, double sigma, double eps);



#endif //MINIMD_HELPER_HPP
