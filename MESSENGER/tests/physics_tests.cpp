#include "../src/trajectory.hpp"
#include "../src/lighttime.hpp"
#include <iostream>

void verify(bool ok, const char* message)
{
    if (!ok)
    {
        throw std::runtime_error(message);
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::string data = argc > 1 ? argv[1] : "data";
        loadKernels(data);
        // Здесь отдельно проверяем базовую ориентацию, без либрации.
        // Рассчитанные амплитуды независимо проверяются libration_tests.
        const auto rotationParameters = loadMercuryParameters(
            std::filesystem::path(data).parent_path() / "config" / "mercury_orbit.txt");
        double zeroLibrations[5]{};
        pdpool_c("BODY199_NUT_PREC_PM", 5, zeroLibrations);
        // Штатный IAU_MERCURY из pck00010.tpc проверяет параметры и формулу.
        // Параметры kernel pool не подменяются числами из проверяемой функции.
        double rotationError = 0;
        for (double t : {0., 354358512., -1e8})
        {
            double own[3][3], spice[3][3];
            mercuryRotation(t, own, rotationParameters, 0.);
            pxform_c("J2000", "IAU_MERCURY", t, spice);
            spiceCheck();
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    rotationError = std::max(rotationError, std::abs(own[i][j] - spice[i][j]));
                }
            }
        }
        verify(rotationError < 2e-12,
               "Mercury Euler rotation differs from independent CSPICE PCK evaluation");
        kclear_c();
        loadKernels(data);
        Gravity gravity(160);
        gravity.load(data + "/HgM009.sha");
        Dynamics dynamics(gravity);
        const double begin = jdTDBtoET(2455646.83);
        const State initial = ephemeris(-236, begin);
        double pericenter = begin, minRadius = norm(pos(initial));
        for (int j = 1; j <= 1200; ++j)
        {
            const double t = begin + 60 * j, r = norm(pos(ephemeris(-236, t)));
            if (r < minRadius)
            {
                minRadius = r;
                pericenter = t;
            }
        }
        const double arcBegin = pericenter - 300;
        const State arcInitial = ephemeris(-236, arcBegin);
        NumericTrajectory dense(dynamics, arcBegin, arcBegin + 600, arcInitial);
        Integrator<Dynamics> direct(dynamics, arcBegin, arcInitial, 5);
        double interpolationError = 0;
        for (int i = 0; i < 20; ++i)
        {
            const double t = arcBegin + 15 + 30 * i;
            interpolationError = std::max(
                interpolationError, 1000 * norm(dense.mercuryRelative(t) - pos(direct.to(t))));
        }
        verify(interpolationError < 1.0,
               "Cubic trajectory interpolation exceeds 1 m on validation arc");
        Stations stations(data);
        SpacecraftSpline reference(begin, begin + 600);
        LightTime light(stations, [&](Epoch t) { return reference(t); });
        const auto path = light.solve(Epoch(begin, .00000001), 34, 34);
        verify(std::max(std::abs(path.up_closure_s), std::abs(path.down_closure_s)) < 1e-8,
               "Light-time equation does not close");
        std::cout << "PASS: independent PCK rotation error " << rotationError
                  << ", cubic interpolation " << interpolationError
                  << " m, split-epoch two-way light time\n";
        kclear_c();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
