#include "Simulation.hpp"

//==================================================================//
// Simulation implementation
//==================================================================//

// Constructor
Simulation::Simulation(unsigned nAtoms, double temp, unsigned runTime) :
    nAtoms_(nAtoms),
    temp_(temp),
    runTime_(runTime)
{
    // Resize the position and velocity matrices to accommodate every atom
    vel_.resize(nAtoms_, 3);
    pos_.resize(nAtoms_, 3);

    // Populate the position matrix with an x, y, and z for each atom
    // by evenly spacing them throughout the cubic box
    // We'll assume a box that is 30 A on a side for now
    int n = (int) (ceil(cbrt(nAtoms_)) + 2);
    double spacing = xVec_[0]/(double (n) - 1.0);
    double xpos = 0.0, ypos = 0.0, zpos = 0.0;

    for (int i = 0, count = 0; i < n-2; i++){
        xpos += spacing;
        for (int j = 0; j < n-2; j++){
            ypos += spacing;
            for (int k = 0; k < n-2; k++){
                zpos += spacing;
                if (count < nAtoms){
                    pos_(count, 0) = xpos;
                    pos_(count, 1) = ypos;
                    pos_(count, 2) = zpos;
                }
                count++;
            }
            zpos = 0.0;
        }
        ypos = 0.0;
    }

    // Populate the velocity matrix with an x, y, and z velocity for each atom
    // Velocities are randomly sampled from a maxwell-boltzmann distribution centered
    // at the appropriate temperature. Therefore individual component velocities
    // are gaussian with mean 0 and std dev sqrt(kT/2m)
    std::default_random_engine generator;
    std::normal_distribution<double> normal(0, sqrt((kb_*temp_)/(2*mass_)));
    for (int i = 0; i < vel_.rows(); i++){
        vel_(i, 0) = normal(generator);
        vel_(i, 1) = normal(generator);
        vel_(i, 2) = normal(generator);
    }
} // Constructor


// Writes current simulation frame to an xyz file
void Simulation::print_pos(){
    long rows = pos_.rows();
    std::ofstream writefile;
    writefile.open ("config.xyz", std::ios::out);
    writefile << rows << "\n#####\n";
    for (int i = 0; i < rows; i++){
        writefile << "Ar\t" << pos_(i, 0) << "\t" << pos_(i, 1) << "\t" << pos_(i, 2) << "\n";
    }
    writefile << std::endl;
    writefile.close();
}

//==================================================================//
// Integrator implementation
//==================================================================//

// Constructor needs a matrix to hold the forces and a second one for accelerations
VelocityVerletIntegrator::VelocityVerletIntegrator(Simulation& sim) : sim_(sim){
    forces.resize(sim_.nAtoms(), 3);
    accel.resize(sim_.nAtoms(), 3);
}


void VelocityVerletIntegrator::Forces(){
    // Parameters
    int nAtoms = sim_.nAtoms();
    double eps = 1.0, sigma = 1.0;
    double dx = 0, dy = 0, dz = 0;
    double dim = sim_.dim();
    double rc = sim_.rc();

    // Zero out the forces, potential, and virial
    for (int i = 0; i < nAtoms; i++){
        forces.row(i) = Eigen::RowVector3d::Zero();
    }
    E_pot = 0;
    virial = 0;

    // Core loop over atom pairs
    for (int i = 0; i < nAtoms; i++){
        for (int j = i+1; j < nAtoms; j++){
            // Distances
            dx = sim_.pos_(i,0) - sim_.pos_(j,0);
            dy = sim_.pos_(i,1) - sim_.pos_(j,1);
            dz = sim_.pos_(i,2) - sim_.pos_(j,2);

            // Minimum image criterion
            if (std::abs(dx) > 0.5*dim) dx -= (dx < 0 ? dim*-1 : dim);
            if (std::abs(dy) > 0.5*dim) dy -= (dy < 0 ? dim*-1 : dim);
            if (std::abs(dz) > 0.5*dim) dz -= (dz < 0 ? dim*-1 : dim);

            // Force from LJ potential
            double len_2 = dx*dx + dy*dy + dz*dz;
            if (len_2 < rc*rc){
                double fr_6 = std::pow(sigma*sigma/len_2, 3);
                double frc = (48.0*eps/len_2)*fr_6*(fr_6 - 0.5);

                // Add individual forces on atoms to force matrix
                // Note force of atom i on atom j is equal to the opposite of the
                // force of j to i by Newton's 3rd law
                forces(i,0) += frc*dx;
                forces(j,0) -= frc*dx;
                forces(i,1) += frc*dy;
                forces(j,1) -= frc*dy;
                forces(i,2) += frc*dz;
                forces(j,2) -= frc*dz;

                // Potential and virial
                E_pot +=  4*eps*fr_6*(fr_6 - 1.0);
                virial += frc*len_2;
            } // if within cutoff
        } // j
    } // i

    // Normalize E_pot and virial
    E_pot /= nAtoms;
    virial /= nAtoms;
} // Forces


void VelocityVerletIntegrator::WrapPos(){
    int nAtoms = sim_.nAtoms();
    double dim = sim_.dim();
    // Loop over the positions of all atoms
    for (int i = 0; i < nAtoms; i++){
        // Loop over x, y, and z directions
        for (int j = 0; j < 3; j++) {
            if (sim_.pos_(i, j) < dim) sim_.pos_(i, j) += dim;
            else if (sim_.pos_(i, j) > dim) sim_.pos_(i, j) -= dim;
        }
    }
}