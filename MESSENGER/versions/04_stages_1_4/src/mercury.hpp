#pragma once
#include "math.hpp"
#include "mercury_parameters.hpp"

inline constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.;

struct MercuryAngles
{
    double rightAscension; // Прямое восхождение полюса alpha, радианы.
    double declination;    // Склонение полюса delta, радианы.
    double primeMeridian;  // Угол начального меридиана W, радианы.
};

// Зависимости alpha(t), delta(t), W(t): Stark et al. (2015), формулы (1), (2).
// Коэффициенты одного набора IAU 2009 перечислены в mercury_parameters.hpp.
inline MercuryAngles mercuryAngles(double et)
{
    const double daysSinceJ2000 = et / 86400.;
    const double centuriesSinceJ2000 = daysSinceJ2000 / 36525.;
    const double alphaDeg = MercuryIAU2009::POLE_RA_DEG +
                           MercuryIAU2009::POLE_RA_RATE_DEG_CENTURY * centuriesSinceJ2000;
    const double deltaDeg = MercuryIAU2009::POLE_DEC_DEG +
                           MercuryIAU2009::POLE_DEC_RATE_DEG_CENTURY * centuriesSinceJ2000;
    double wDeg = MercuryIAU2009::PRIME_MERIDIAN_DEG +
                  MercuryIAU2009::SPIN_RATE_DEG_DAY * daysSinceJ2000;
    for (const auto& term : MercuryIAU2009::LIBRATIONS)
    {
        const double argument = (term.phaseDeg + term.rateDegDay * daysSinceJ2000) * DEG_TO_RAD;
        wDeg += term.amplitudeDeg * std::sin(argument);
    }
    return {alphaDeg * DEG_TO_RAD, deltaDeg * DEG_TO_RAD,
            std::remainder(wDeg, 360.) * DEG_TO_RAD};
}

// Произведение матриц 3x3; временная матрица допускает совпадение входа и выхода.
inline void multiplyRotations(const double left[3][3], const double right[3][3],
                              double result[3][3])
{
    double product[3][3]{};
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            for (int k = 0; k < 3; ++k)
            {
                product[row][col] += left[row][k] * right[k][col];
            }
        }
    }
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            result[row][col] = product[row][col];
        }
    }
}

// Небесные оси -> оси Меркурия: композиция трех поворотов, Stark (2015), (5).
// Для обычных пассивных R1/R3 и астрономических alpha/delta:
// R = R3(W) * R1(pi/2 - delta) * R3(pi/2 + alpha).
// Прямое восхождение и склонение полюса не являются самими углами Эйлера.
inline void mercuryRotation(double et, double rotation[3][3])
{
    const auto angles = mercuryAngles(et);
    const double poleAzimuth = 90 * DEG_TO_RAD + angles.rightAscension;
    const double poleTilt = 90 * DEG_TO_RAD - angles.declination;
    const double spin = angles.primeMeridian;

    const double ca = std::cos(poleAzimuth), sa = std::sin(poleAzimuth);
    const double cd = std::cos(poleTilt), sd = std::sin(poleTilt);
    const double cw = std::cos(spin), sw = std::sin(spin);

    // Первый поворот вокруг Z задаёт азимут направления полюса.
    const double aroundZ[3][3] = {{ca, sa, 0}, {-sa, ca, 0}, {0, 0, 1}};
    // Второй поворот вокруг X совмещает ось Z с полюсом Меркурия.
    const double aroundX[3][3] = {{1, 0, 0}, {0, cd, sd}, {0, -sd, cd}};
    // Третий поворот вокруг новой оси Z учитывает собственное вращение.
    const double aroundPole[3][3] = {{cw, sw, 0}, {-sw, cw, 0}, {0, 0, 1}};

    double poleOrientation[3][3];
    multiplyRotations(aroundX, aroundZ, poleOrientation);
    multiplyRotations(aroundPole, poleOrientation, rotation);
}
