#pragma once
#include "spice.hpp"
#include "integrator.hpp"
#include <vector>

class NumericTrajectory
{
    struct Node
    {
        double t;
        State state;
        Vec acceleration;
    };

    std::vector<Node> nodes_;

public:
    NumericTrajectory(const Dynamics& dynamics, double begin, double end, State initial,
                      double step = 60)
    {
        if (end <= begin)
        {
            throw std::runtime_error("Empty trajectory interval");
        }
        Integrator<Dynamics> integrator(dynamics, begin, initial, step);
        for (double t = begin; t <= end + step; t += step)
        {
            const auto s = integrator.to(t);
            nodes_.push_back({t, s, dynamics.acceleration(t, s)});
        }
    }

    Vec mercuryRelative(const Epoch& et) const
    {
        const double epoch = et.high;
        if (epoch < nodes_.front().t || epoch > nodes_.back().t)
        {
            throw std::runtime_error("Trajectory query outside integrated interval");
        }
        auto found = std::upper_bound(nodes_.begin(), nodes_.end(), epoch,
                                      [](double t, const Node& n) { return t < n.t; });
        size_t i = found == nodes_.end() ? nodes_.size() - 2 : size_t(found - nodes_.begin() - 1);
        const auto& a = nodes_[i];
        const auto& b = nodes_[i + 1];
        const double h = b.t - a.t, u = ((et.high - a.t) + et.low) / h;
        Vec r{};
        for (int j = 0; j < 3; ++j)
        {
            const double v0 = h * a.state[j + 3], a0 = h * h * a.acceleration[j];
            const double A = b.state[j] - a.state[j] - v0 - .5 * a0,
                         B = h * b.state[j + 3] - v0 - a0, C = h * h * b.acceleration[j] - a0;
            r[j] = a.state[j] + u * (v0 + u * (.5 * a0 + u * ((10 * A - 4 * B + .5 * C) +
                                                              u * ((-15 * A + 7 * B - C) +
                                                                   u * (6 * A - 3 * B + .5 * C)))));
        }
        return r;
    }

    Vec operator()(const Epoch& et) const
    {
        return pos(ephemeris(199, et, 0)) + mercuryRelative(et);
    }
};
