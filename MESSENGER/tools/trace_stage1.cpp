#include "spice.hpp"
#include "integrator.hpp"
int main(int argc, char** argv)
{
    loadKernels(argv[1]);
    Gravity gravity(160);
    gravity.load(std::string(argv[1]) + "/HgM009.sha");
    Dynamics dynamics(gravity);
    const double t0 = jdTDBtoET(2455646.83);
    Integrator<Dynamics> integrator(dynamics, t0, ephemeris(-236, t0), 60);
    integrator.to(t0 + 600);
}
