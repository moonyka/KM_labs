#pragma once
#include "math.hpp"
#include <vector>
inline constexpr double LIGHT_SPEED_KM_S = 299792.458;

struct GravBody
{
    int id;
    double gm;
    Vec r, v, a{};
    double potential = 0;
};

// Ньютоновские ускорения и потенциалы массивных тел для релятивистской поправки.
inline void prepareEIH(std::vector<GravBody>& bodies)
{
    for (size_t j = 0; j < bodies.size(); ++j)
    {
        auto& b = bodies[j];
        b.a = {};
        b.potential = 0;
        for (size_t k = 0; k < bodies.size(); ++k)
        {
            if (k != j)
            {
                const Vec d = bodies[k].r - b.r;
                const double r = norm(d);
                b.potential += bodies[k].gm / r;
                b.a = b.a + (bodies[k].gm / (r * r * r)) * d;
            }
        }
    }
}

// Поправка первого постньютоновского приближения в общей барицентрической системе.
inline Vec eihCorrection(const Vec& r, const Vec& v, const std::vector<GravBody>& bodies,
                         int excluded = -1)
{
    double potential = 0;
    for (const auto& b : bodies)
    {
        if (b.id != excluded)
        {
            potential += b.gm / norm(b.r - r);
        }
    }
    Vec correction{};
    const double invC2 = 1 / (LIGHT_SPEED_KM_S * LIGHT_SPEED_KM_S);
    for (const auto& b : bodies)
    {
        if (b.id != excluded)
        {
            const Vec d = b.r - r;
            const double distance = norm(d), invR3 = 1 / (distance * distance * distance),
                         projected = dot(d, b.v) / distance;
            const double scalar = dot(v, v) + 2 * dot(b.v, b.v) - 4 * dot(v, b.v) -
                                  1.5 * projected * projected - 4 * potential - b.potential +
                                  .5 * dot(d, b.a);
            correction = correction + (b.gm * invR3 * scalar) * d;
            correction = correction + (-b.gm * invR3 * dot(d, 4. * v - 3. * b.v)) * (v - b.v);
            correction = correction + (3.5 * b.gm / distance) * b.a;
        }
    }
    return invC2 * correction;
}

