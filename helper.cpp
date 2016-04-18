//
//

#include "helper.hpp"


Eigen::RowVector3d LJForce(Eigen::RowVector3d dist_v, double rc, double sigma = 3.345, double eps = 1.654e-21){
    Eigen::RowVector3d force_v;
    double dist_norm = dist_v.norm();
    if (dist_norm > rc){
        return Eigen::RowVector3d::Zero().eval(); // Return 0 force vector if distance is greater than cutoff radius
    }
    force_v = (dist_v/dist_norm) * (24*eps/sigma) * (2*pow(sigma/dist_norm, 13) - pow(sigma/dist_norm, 7));
    return force_v;
}


Eigen::RowVector3d PeriodicDistanceVec(Eigen::RowVector3d x_pos, Eigen::RowVector3d y_pos, double dim){
    Eigen::RowVector3d dist_v;
    for (int i = 0; i < 3; i++){
        double diff = x_pos(i)-y_pos(i);
        if (diff > 0.5*dim) dist_v(i) = dim - diff; // Account for periodic boundary conditions
        else dist_v(i) = diff;
    }
    return dist_v;
}

