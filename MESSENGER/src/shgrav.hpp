#pragma once
#include "gravity.hpp"
inline constexpr double AU_KM = 149597870.7, DAY_SECONDS = 86400.;
#ifdef MESSENGER_USE_SHTOOLS
#include "grav.h"
#else
struct ShGrav
{
    int lmax;
    double gm, r0;
    Gravity field;
};

inline ShGrav shGravCreate(int lmax, double gm, double radius)
{
    ShGrav value{lmax, gm, radius, Gravity(lmax)};
    value.field.gm = gm * AU_KM * AU_KM * AU_KM / (DAY_SECONDS * DAY_SECONDS);
    value.field.radius = radius * AU_KM;
    return value;
}

inline int shGravLoadCoefficients(ShGrav& field, const char* path)
{
    field.field.load(path);
    return 1;
}

inline void shGravCalculate(const ShGrav& field, double x, double y, double z, double result[3])
{
    const Vec r = AU_KM * Vec{x, y, z};
    const double distance = norm(r);
    const Vec central = (-field.field.gm / (distance * distance * distance)) * r;
    const Vec harmonics = field.field.acceleration(r) - central;
    for (int i = 0; i < 3; ++i)
    {
        result[i] = harmonics[i] * DAY_SECONDS * DAY_SECONDS / AU_KM;
    }
}
#endif
