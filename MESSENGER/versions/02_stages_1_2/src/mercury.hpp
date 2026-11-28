#pragma once
#include "math.hpp"
#include "libration.hpp"

inline constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.;

struct MercuryAngles
{
    double rightAscension; // Прямое восхождение полюса alpha, радианы.
    double declination;    // Склонение полюса delta, радианы.
    double primeMeridian;  // Угол начального меридиана W, радианы.
};

// Физические исходные параметры читаются из PCK и файла орбитальных данных.
// relativeAmplitudes содержит результат формулы (3) Stark, а не готовую таблицу.
struct MercuryParameters
{
    std::array<double, 3> poleRA{}, poleDec{}, primeMeridian{};
    double phaseAtJ2000Deg = 0;
    double meanMotionDegDay = 0;
    double fundamentalLibrationDeg = 0;
    double eccentricity = 0;
    std::array<double, 5> relativeAmplitudes{};
};

// Stark et al. (2015), (1), (2): углы вычисляются для текущей эпохи.
// В PCK полюс задан полиномами по векам, начальный меридиан — по суткам.
inline MercuryAngles mercuryAngles(double et, const MercuryParameters& model,
                                   double librationScale = 1.)
{
    const double daysSinceJ2000 = et / 86400.;
    const double centuriesSinceJ2000 = daysSinceJ2000 / 36525.;
    const auto polynomial = [](const std::array<double, 3>& a, double t)
    {
        return a[0] + t * (a[1] + t * a[2]);
    };
    const double alphaDeg = polynomial(model.poleRA, centuriesSinceJ2000);
    const double deltaDeg = polynomial(model.poleDec, centuriesSinceJ2000);
    double wDeg = polynomial(model.primeMeridian, daysSinceJ2000);
    // Один исходный аргумент задаёт все k гармоник: k*n0*(d+t0).
    for (int k = 1; k <= 5; ++k)
    {
        const double argument = k * (model.phaseAtJ2000Deg +
                                     model.meanMotionDegDay * daysSinceJ2000) * DEG_TO_RAD;
        wDeg += librationScale * model.fundamentalLibrationDeg *
                model.relativeAmplitudes[k - 1] * std::sin(argument);
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
inline void mercuryRotation(double et, double rotation[3][3],
                            const MercuryParameters& model, double librationScale = 1.)
{
    const auto angles = mercuryAngles(et, model, librationScale);
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
