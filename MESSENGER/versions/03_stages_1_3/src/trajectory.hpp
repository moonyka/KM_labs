#pragma once
#include "spice.hpp"
#include "integrator.hpp"
#include "spline.hpp"
#include <vector>

// Численная орбита первого этапа; между узлами — натуральный кубический сплайн.
// Burden & Faires, Numerical Analysis, 9 ed., раздел 3.5, алгоритм 3.4, стр. 149.
class NumericTrajectory
{
    double begin_, end_;
    std::array<CubicSpline, 3> coordinates_;

public:
    NumericTrajectory(const Dynamics& dynamics, double begin, double end, State initial,
                      double step = 60)
        : begin_(begin), end_(end)
    {
        if (end <= begin)
        {
            throw std::runtime_error("Empty trajectory interval");
        }
        const int paddingSteps = 32;
        std::vector<double> times;
        std::array<std::vector<double>, 3> positions;
        const auto append = [&](double t, const State& state)
        {
            times.push_back(t - begin);
            for (int j = 0; j < 3; ++j)
            {
                positions[j].push_back(state[j]);
            }
        };

        // Дополнительные узлы уносят искусственные условия S''=0 за рабочий интервал.
        // В обратном направлении также используется тот же RK8, без эталонных координат.
        Integrator<Dynamics> backward(dynamics, begin, initial, step);
        std::vector<State> before;
        for (int i = 1; i <= paddingSteps; ++i)
        {
            before.push_back(backward.to(begin - i * step));
        }
        for (int i = paddingSteps; i > 0; --i)
        {
            append(begin - i * step, before[i - 1]);
        }
        // Основное интегрирование стартует с точного переданного начального состояния.
        Integrator<Dynamics> forward(dynamics, begin, initial, step);
        append(begin, initial);
        for (double t = begin + step; t <= end + paddingSteps * step; t += step)
        {
            append(t, forward.to(t));
        }
        for (int j = 0; j < 3; ++j)
        {
            coordinates_[j] = CubicSpline(times, positions[j]);
        }
    }

    Vec mercuryRelative(const Epoch& et) const
    {
        const double t = (et.high - begin_) + et.low;
        if (t < 0 || t > end_ - begin_)
        {
            throw std::runtime_error("Trajectory query outside integrated interval");
        }
        return {coordinates_[0](t), coordinates_[1](t), coordinates_[2](t)};
    }

    Vec operator()(const Epoch& et) const
    {
        return pos(ephemeris(199, et, 0)) + mercuryRelative(et);
    }
};
