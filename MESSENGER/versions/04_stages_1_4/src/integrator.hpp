#pragma once
#include "math.hpp"
#include "rk8_table.hpp"

struct IntegrationStats
{
    size_t accepted = 0;
    size_t rejected = 0;
};

template <class F>
class Integrator
{
    const F& derivative_;
    double time_;
    double step_;
    State state_;

public:
    IntegrationStats stats;

    Integrator(const F& derivative, double initialTime, State initialState, double step = 60)
        : derivative_(derivative), time_(initialTime), step_(step), state_(initialState)
    {
        if (!std::isfinite(step) || step <= 0)
        {
            throw std::runtime_error("Invalid RK8 step");
        }
    }

    const State& state() const
    {
        return state_;
    }

    State to(double target)
    {
        const double direction = target >= time_ ? 1. : -1.;

        while (direction * (target - time_) > 1e-7)
        {
            const double h = direction * std::min(step_, std::abs(target - time_));
            State stageDerivatives[RK_STAGES]{};

            // k_i = f(t + c_i*h, s + h*sum(a_ij*k_j)).
            for (int stage = 0; stage < RK_STAGES; ++stage)
            {
                State stageState = state_;

                for (int component = 0; component < 6; ++component)
                {
                    double derivativeSum = 0;
                    for (int previousStage = 0; previousStage < stage; ++previousStage)
                    {
                        derivativeSum +=
                            RK_A[stage][previousStage] * stageDerivatives[previousStage][component];
                    }
                    stageState[component] += h * derivativeSum;
                }

                stageDerivatives[stage] = derivative_(time_ + RK_C[stage] * h, stageState);
            }

            // s_next = s + h*sum(b_i*k_i).
            for (int component = 0; component < 6; ++component)
            {
                double derivativeSum = 0;
                for (int stage = 0; stage < RK_STAGES; ++stage)
                {
                    derivativeSum += RK_B[stage] * stageDerivatives[stage][component];
                }
                state_[component] += h * derivativeSum;

                if (!std::isfinite(state_[component]))
                {
                    throw std::runtime_error("Nonfinite RK8 state");
                }
            }

            time_ += h;
            ++stats.accepted;
        }

        time_ = target;
        return state_;
    }
};
