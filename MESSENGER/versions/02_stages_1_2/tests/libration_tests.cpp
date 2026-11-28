#include "../src/libration.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

void verifyLibration(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

int main(int argc, char** argv)
{
    try
    {
        verifyLibration(argc == 2, "Missing independent libration reference");
        std::ifstream input(argv[1]);
        verifyLibration(bool(input), "Cannot open libration reference");
        std::string row;
        std::getline(input, row);
        int checked = 0;
        double maxError = 0.;
        while (std::getline(input, row))
        {
            std::replace(row.begin(), row.end(), ',', ' ');
            std::istringstream fields(row);
            double eccentricity, expected;
            int q;
            fields >> eccentricity >> q >> expected;
            verifyLibration(bool(fields), "Invalid independent reference row");
            const double error = std::abs(kaulaG20(q, eccentricity) - expected);
            maxError = std::max(maxError, error);
            verifyLibration(error < 2e-13, "Kaula integral differs from independent mean-anomaly quadrature");
            ++checked;
        }
        verifyLibration(checked == 44, "Incomplete independent libration reference");
        const auto circular = calculateLibrationRatios(0.);
        verifyLibration(std::abs(circular[0]-1) < 1e-15, "Wrong circular fundamental");
        for (int k = 1; k < 5; ++k)
            verifyLibration(std::abs(circular[k]) < 1e-14, "Spurious higher circular harmonics");
        const auto coarse = calculateLibrationRatios(.2056317, 1024);
        const auto fine = calculateLibrationRatios(.2056317, 4096);
        for (int k = 0; k < 5; ++k)
            verifyLibration(std::abs(coarse[k] - fine[k]) < 2e-13, "Libration quadrature has not converged");
        const auto other = calculateLibrationRatios(.1);
        verifyLibration(std::abs(other[1] - fine[1]) > .01, "Libration coefficients do not depend on eccentricity");
        bool invalid = false;
        try { kaulaG20(0, 1.); } catch (const std::runtime_error&) { invalid = true; }
        verifyLibration(invalid, "Invalid eccentricity accepted");
        std::cout << "PASS: 44 independent Kaula integrals, circular limit, refinement, changing eccentricity; max error " << maxError << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
