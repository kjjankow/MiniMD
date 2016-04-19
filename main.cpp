#include "Simulation.hpp"

int main() {
    Simulation ArgonSim = Simulation(125, 300, 1000);
    VelocityVerletIntegrator integrator(ArgonSim);
    integrator.Forces();
    for (int i = 0; i < 101; i++) {
        ArgonSim.incTime();
        integrator.VelocityVerletStep(0.0005);
        integrator.SetValues();
        ArgonSim.print_stat();
        if (i%10==0) ArgonSim.print_pos();
    }

    std::cout << "The simulation has completed successfully!" << std::endl;
    std::cout << "Results have been written to config.xyz" << std::endl;
    return 0;
}