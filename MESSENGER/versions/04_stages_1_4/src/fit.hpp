#pragma once
#include "doppler.hpp"
#include "trajectory.hpp"

inline State solveNormal(std::array<State, 6> matrix, State rhs)
{
    for (int col = 0; col < 6; ++col)
    {
        int pivot = col;
        for (int row = col + 1; row < 6; ++row)
        {
            if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col]))
            {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][col]) < 1e-20)
        {
            throw std::runtime_error("Singular least-squares system");
        }
        std::swap(matrix[pivot], matrix[col]);
        std::swap(rhs[pivot], rhs[col]);
        double d = matrix[col][col];
        for (int j = col; j < 6; ++j)
        {
            matrix[col][j] /= d;
        }
        rhs[col] /= d;
        for (int row = 0; row < 6; ++row)
        {
            if (row != col)
            {
                d = matrix[row][col];
                for (int j = col; j < 6; ++j)
                {
                    matrix[row][j] -= d * matrix[col][j];
                }
                rhs[row] -= d * rhs[col];
            }
        }
    }
    return rhs;
}

inline double rms(const std::vector<double>& values)
{
    double sum = 0;
    for (double x : values)
    {
        sum += x * x;
    }
    return std::sqrt(sum / values.size());
}

inline void runFit(const std::string& data, const std::string& out, int points, int iterations,
                   double fitHours, double countScale = .001, double integrationStep = 60)
{
    auto all = readObservations(data + "/obs_4.txt", countScale);
    all.erase(std::remove_if(all.begin(), all.end(), [](const auto& o) { return !o.eligible; }),
              all.end());
    if (fitHours > 0)
    {
        const double end = all.front().tai + fitHours * 3600;
        all.erase(
            std::remove_if(all.begin(), all.end(), [&](const auto& o) { return o.tai > end; }),
            all.end());
    }
    const size_t stride = std::max(size_t(1), (all.size() + points - 1) / points);
    std::vector<Observation> obs;
    for (size_t i = 0; i < all.size(); i += stride)
    {
        obs.push_back(all[i]);
    }
    if (obs.size() < 12)
    {
        throw std::runtime_error("At least 12 fit observations required");
    }
    const Stations stations(data);
    const Ramps ramps(data + "/ramps_4.txt");
    Gravity gravity(160);
    gravity.load(data + "/HgM009.sha");
    Dynamics dynamics(gravity);
    const double begin = taiOffsetToET(obs.front().tai) - 1800,
                 end = taiOffsetToET(obs.back().tai) + 120;
    const State reference = ephemeris(-236, begin);
    State parameters{};
    const State scales{1, 1, 1, .001, .001, .001};
    const auto evaluate = [&](const State& p)
    {
        State initial = reference;
        for (int j = 0; j < 6; ++j)
        {
            initial[j] += p[j] * scales[j];
        }
        NumericTrajectory trajectory(dynamics, begin, end, initial, integrationStep);
        LightTime light(stations, [&](Epoch et) { return trajectory(et); });
        std::vector<double> result;
        for (const auto& o : obs)
        {
            result.push_back(computeDoppler(o, ramps, light, stations).computed);
        }
        return result;
    };
    const auto residuals = [&](const std::vector<double>& computed)
    {
        std::vector<double> result;
        for (size_t i = 0; i < obs.size(); ++i)
        {
            result.push_back(obs[i].observed - computed[i]);
        }
        return result;
    };
    std::filesystem::create_directories(out);
    std::ofstream history(out + "/iterations.csv");
    history << "iteration,attempt,rms_hz,accepted,lambda\n";
    std::cout << "Fitting " << obs.size() << " observations over " << (end - begin) / 3600
              << " h; initial orbit integration..." << std::endl;
    auto computed = evaluate(parameters);
    auto residual = residuals(computed);
    const double initialRms = rms(residual);
    double error = initialRms, lambda = .001;
    history << std::setprecision(16) << 0 << ",0," << error << ",1," << lambda << '\n';
    history.flush();
    std::cout << "Initial RMS " << error << " Hz" << std::endl;
    bool converged = false;
    int accepted = 0;
    for (int iteration = 1; iteration <= iterations; ++iteration)
    {
        std::vector<State> jacobian(obs.size());
        for (int j = 0; j < 6; ++j)
        {
            std::cout << "Iteration " << iteration << ": derivative " << j + 1 << "/6" << std::endl;
            State shifted = parameters;
            shifted[j] += .01;
            const auto values = evaluate(shifted);
            for (size_t i = 0; i < obs.size(); ++i)
            {
                jacobian[i][j] = (values[i] - computed[i]) / .01;
            }
        }
        std::array<State, 6> normal{};
        State rhs{};
        for (size_t i = 0; i < obs.size(); ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                rhs[j] += jacobian[i][j] * residual[i];
                for (int k = 0; k < 6; ++k)
                {
                    normal[j][k] += jacobian[i][j] * jacobian[i][k];
                }
            }
        }
        bool improved = false;
        double improvement = 0, stepNorm = 0;
        for (int attempt = 0; attempt < 8; ++attempt)
        {
            auto damped = normal;
            for (int j = 0; j < 6; ++j)
            {
                damped[j][j] += lambda * std::max(normal[j][j], 1e-12);
            }
            State step = solveNormal(damped, rhs), trial = parameters;
            stepNorm = 0;
            for (int j = 0; j < 6; ++j)
            {
                step[j] = std::clamp(step[j], -10., 10.);
                trial[j] += step[j];
                stepNorm = std::max(stepNorm, std::abs(step[j]));
            }
            std::cout << "Iteration " << iteration << ": trial, lambda=" << lambda << std::endl;
            auto trialComputed = evaluate(trial);
            auto trialResidual = residuals(trialComputed);
            const double trialRms = rms(trialResidual);
            history << iteration << ',' << attempt + 1 << ',' << trialRms << ','
                    << (trialRms < error) << ',' << lambda << '\n';
            history.flush();
            if (trialRms < error)
            {
                improvement = error - trialRms;
                error = trialRms;
                parameters = trial;
                computed = std::move(trialComputed);
                residual = std::move(trialResidual);
                lambda = std::max(lambda / 3, 1e-9);
                improved = true;
                ++accepted;
                break;
            }
            lambda *= 10;
        }
        std::cout << "Iteration " << iteration << ": RMS " << error << " Hz, accepted=" << improved
                  << std::endl;
        if (improved && (stepNorm < 1e-5 || improvement < 1e-5))
        {
            converged = true;
            break;
        }
        if (!improved)
        {
            break;
        }
    }
    std::ofstream f(out + "/fit.csv");
    f << std::setprecision(16) << "tai_offset_s,station,observed_hz,computed_hz,residual_hz\n";
    for (size_t i = 0; i < obs.size(); ++i)
    {
        f << obs[i].tai << ',' << obs[i].rx << ',' << obs[i].observed << ',' << computed[i] << ','
          << residual[i] << '\n';
    }
    std::ofstream summary(out + "/summary.json");
    summary << std::setprecision(16)
            << "{\n  \"method\": \"Levenberg-Marquardt, unweighted least squares\",\n  "
               "\"parameters\": \"six initial Cartesian state corrections\",\n  \"epoch_et\": "
            << begin << ",\n  \"observation_count\": " << obs.size()
            << ",\n  \"fit_hours_option\": " << fitHours
            << ",\n  \"count_scale_seconds\": " << countScale
            << ",\n  \"integration_step_s\": " << integrationStep
            << ",\n  \"trajectory_interpolation\": \"quintic Hermite\",\n  \"initial_rms_hz\": "
            << initialRms << ",\n  \"final_rms_hz\": " << error
            << ",\n  \"accepted_iterations\": " << accepted
            << ",\n  \"converged\": " << (converged ? "true" : "false")
            << ",\n  \"correction_km_and_km_s\": [";
    for (int j = 0; j < 6; ++j)
    {
        if (j)
        {
            summary << ", ";
        }
        summary << parameters[j] * scales[j];
    }
    summary << "],\n  \"fitted_state_km_and_km_s\": [";
    for (int j = 0; j < 6; ++j)
    {
        if (j)
        {
            summary << ", ";
        }
        summary << reference[j] + parameters[j] * scales[j];
    }
    summary << "]\n}\n";
    std::cout << "Fit complete: RMS " << initialRms << " -> " << error
              << " Hz; converged=" << converged << std::endl;
}
