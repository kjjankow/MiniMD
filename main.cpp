#include "Simulation.hpp"

int main() {
    Simulation ArgonSim = Simulation(125, 300, 1000);
    VelocityVerletIntegrator integrator(ArgonSim);
    integrator.Forces();
    integrator.WrapPos();
    ArgonSim.print_pos();
    std::cout << "The simulation has completed successfully!" << std::endl;
    std::cout << "Results have been written to config.xyz" << std::endl;
    return 0;
}