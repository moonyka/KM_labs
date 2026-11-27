#pragma once
#include <array>

// Единый набор IAU 2009: Archinal et al. (2011), таблица 1, строка Mercury.
// DOI: 10.1007/s10569-010-9320-4 (первая статья из письма от 23.10).
// Числа также опубликованы NAIF в pck00010.tpc, раздел Mercury.
namespace MercuryIAU2009
{
inline constexpr char NAME[] = "IAU 2009 / Archinal 2011";
inline constexpr double POLE_RA_DEG = 281.0097;
inline constexpr double POLE_RA_RATE_DEG_CENTURY = -0.0328;
inline constexpr double POLE_DEC_DEG = 61.4143;
inline constexpr double POLE_DEC_RATE_DEG_CENTURY = -0.0049;
inline constexpr double PRIME_MERIDIAN_DEG = 329.5469;
inline constexpr double SPIN_RATE_DEG_DAY = 6.1385025;

struct LibrationTerm
{
    double amplitudeDeg; // Амплитуда слагаемого угла W, градусы.
    double phaseDeg;     // Фаза аргумента синуса на J2000, градусы.
    double rateDegDay;   // Изменение аргумента за сутки, градусы/сутки.
};

// W = W0 + W1*d + sum(Ak*sin(Mk0 + Mk1*d)).
// Опубликованные амплитуды и фазы, без перенормировки на 40 секунд дуги.
inline constexpr std::array<LibrationTerm, 5> LIBRATIONS{{
    {0.00993822, 174.791086, 4.092335},
    {-0.00104581, 349.582171, 8.184670},
    {-0.00010280, 164.373257, 12.277005},
    {-0.00002364, 339.164343, 16.369340},
    {-0.00000532, 153.955429, 20.461675},
}};
}
