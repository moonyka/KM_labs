#pragma once
#include "math.hpp"
#include <vector>

inline constexpr double LIGHT_SPEED_KM_S = 299792.458;
inline constexpr double PPN_BETA = 1.;
inline constexpr double PPN_GAMMA = 1.;

struct GravBody
{
    int id;
    double gm;
    Vec r, v;                      // Барицентрические положение и скорость: км, км/с.
    Vec newtonianAcceleration{};   // Ускорение тела j от остальных массивных тел: км/с².
    double externalPotential = 0.; // U_j = sum(k != j) GM_k/r_jk: км²/с², положительная сумма.
};

// Мойер, раздел 4.4.1, (4-27), стр. 4-20; пояснение к (4-26), стр. 4-21.
// Для каждого массивного тела j вычисляем a_j^N и U_j в текущую эпоху.
// Они нужны внутри поправок 1/c²; это не повторное добавление сил к аппарату.
// Вызов: Dynamics::environmentAt -> prepareEIH -> eihCorrection для Меркурия и КА.
inline void prepareEIH(std::vector<GravBody>& bodies)
{
    for (size_t j = 0; j < bodies.size(); ++j)
    {
        auto& bodyJ = bodies[j];
        bodyJ.newtonianAcceleration = {};
        bodyJ.externalPotential = 0.;
        for (size_t k = 0; k < bodies.size(); ++k)
        {
            if (k == j)
            {
                continue; // Собственное поле в точке массы не вычисляется.
            }
            const Vec separationKJ = bodies[k].r - bodyJ.r;
            const double distanceKJ = norm(separationKJ);
            if (distanceKJ == 0.)
            {
                throw std::runtime_error("Coincident massive bodies in EIH preparation");
            }
            bodyJ.externalPotential += bodies[k].gm / distanceKJ;
            bodyJ.newtonianAcceleration = bodyJ.newtonianAcceleration +
                (bodies[k].gm / (distanceKJ * distanceKJ * distanceKJ)) * separationKJ;
        }
    }
}

// Релятивистская часть (4-26) Мойера, стр. 4-19, при beta = gamma = 1.
// Ньютоновский член сюда не входит: он отдельно вычисляется в Dynamics::acceleration.
// excludedId исключает собственную массу, когда целевое тело i — Меркурий.
inline Vec eihCorrection(const Vec& positionI, const Vec& velocityI,
                         const std::vector<GravBody>& bodies, int excludedId = -1)
{
    double potentialI = 0.;
    for (const auto& bodyJ : bodies)
    {
        if (bodyJ.id != excludedId)
        {
            potentialI += bodyJ.gm / norm(bodyJ.r - positionI);
        }
    }

    Vec correction{};
    for (const auto& bodyJ : bodies)
    {
        if (bodyJ.id == excludedId)
        {
            continue;
        }
        const Vec separationJI = bodyJ.r - positionI;
        const double distanceIJ = norm(separationJI);
        const double gmOverR3 = bodyJ.gm / (distanceIJ * distanceIJ * distanceIJ);
        const double projectedVelocityJ = dot(separationJI, bodyJ.v) / distanceIJ;

        // Скаляр в фигурных скобках первой суммы (4-26), без единицы Ньютона.
        const double potentialTerms = -2 * (PPN_BETA + PPN_GAMMA) * potentialI -
                                      (2 * PPN_BETA - 1) * bodyJ.externalPotential;
        const double velocityTerms = PPN_GAMMA * dot(velocityI, velocityI) +
            (1 + PPN_GAMMA) * dot(bodyJ.v, bodyJ.v) -
            2 * (1 + PPN_GAMMA) * dot(velocityI, bodyJ.v) -
            1.5 * projectedVelocityJ * projectedVelocityJ;
        const double accelerationTerm = 0.5 * dot(separationJI, bodyJ.newtonianAcceleration);
        correction = correction +
            (gmOverR3 * (potentialTerms + velocityTerms + accelerationTerm)) * separationJI;

        // Вторая сумма (4-26); минус возникает из r_i-r_j = -separationJI.
        const Vec weightedVelocity = (2 + 2 * PPN_GAMMA) * velocityI -
                                     (1 + 2 * PPN_GAMMA) * bodyJ.v;
        correction = correction - (gmOverR3 * dot(separationJI, weightedVelocity)) *
                                  (velocityI - bodyJ.v);

        // Последняя сумма (4-26): (3+4*gamma)/2 * sum(GM_j*a_j/r_ij).
        correction = correction + ((3 + 4 * PPN_GAMMA) * 0.5 * bodyJ.gm / distanceIJ) *
                                  bodyJ.newtonianAcceleration;
    }
    return (1 / (LIGHT_SPEED_KM_S * LIGHT_SPEED_KM_S)) * correction;
}
