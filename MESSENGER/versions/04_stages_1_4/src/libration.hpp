#pragma once
#include <array>
#include <cmath>
#include <stdexcept>

// G_{20q}(e) = X^{-3,2}_{2+q}(e), функция эксцентриситета Каулы.
// Из определения коэффициента Фурье Хансена переходим к эксцентрической аномалии E:
// G = (1/2pi) * integral[0,2pi] cos(2f(E)-(2+q)*M(E))/(1-e*cos(E))^2 dE.
// M=E-e*sin(E), cos(f)=(cos(E)-e)/(1-e*cos(E)), sin(f)=sqrt(1-e^2)*sin(E)/(1-e*cos(E)).
// Интеграл вычисляется составной формулой Симпсона (Numerical Analysis, раздел 4.4).
inline double kaulaG20(int q, double eccentricity, int subdivisions = 4096)
{
    if (!std::isfinite(eccentricity) || eccentricity < 0 || eccentricity >= 1 ||
        subdivisions < 2 || subdivisions % 2 != 0)
    {
        throw std::runtime_error("Invalid Kaula quadrature parameters");
    }
    constexpr double pi = 3.14159265358979323846;
    const double h = 2 * pi / subdivisions;
    const double minorAxisRatio = std::sqrt(1 - eccentricity * eccentricity);
    double sum = 0., compensation = 0.;
    for (int j = 0; j <= subdivisions; ++j)
    {
        const double E = j * h;
        const double denominator = 1 - eccentricity * std::cos(E);
        const double cosF = (std::cos(E) - eccentricity) / denominator;
        const double sinF = minorAxisRatio * std::sin(E) / denominator;
        const double M = E - eccentricity * std::sin(E);
        const double phase = (2 + q) * M;
        const double integrand = ((cosF * cosF - sinF * sinF) * std::cos(phase) +
                                  2 * cosF * sinF * std::sin(phase)) /
                                 (denominator * denominator);
        const int weight = (j == 0 || j == subdivisions) ? 1 : (j % 2 ? 4 : 2);
        // Компенсированное суммирование уменьшает потерю точности от взаимного вычитания.
        const double term = weight * integrand - compensation;
        const double next = sum + term;
        compensation = (next - sum) - term;
        sum = next;
    }
    return (h / (6 * pi)) * sum;
}

// Stark et al. (2015), (3): F_k=(G_{20,1-k}-G_{20,1+k})/k^2;
// g_(88/(k+1)) = g_(88/k) * F_(k+1)/F_k.
// Для численной устойчивости та же рекурсия записана как g_k/g_1 = F_k/F_1.
inline std::array<double, 5> calculateLibrationRatios(double eccentricity,
                                                    int subdivisions = 4096)
{
    std::array<double, 5> factors{};
    for (int k = 1; k <= 5; ++k)
    {
        factors[k - 1] = (kaulaG20(1 - k, eccentricity, subdivisions) -
                          kaulaG20(1 + k, eccentricity, subdivisions)) / (k * k);
    }
    if (std::abs(factors[0]) < 1e-12)
    {
        throw std::runtime_error("Vanishing fundamental libration factor");
    }
    const double fundamental = factors[0];
    for (auto& factor : factors)
    {
        factor /= fundamental;
    }
    return factors;
}
