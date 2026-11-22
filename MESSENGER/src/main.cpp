#include "gravity.hpp"
#include "spice.hpp"
#include "integrator.hpp"
#ifndef MESSENGER_STAGE
#define MESSENGER_STAGE 4
#endif
#if MESSENGER_STAGE >= 2
#include "lighttime.hpp"
#include "trajectory.hpp"
#endif
#if MESSENGER_STAGE >= 3
#include "doppler.hpp"
#endif
#if MESSENGER_STAGE >= 4
#include "fit.hpp"
#endif
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <optional>

void require(bool ok, const char* message)
{
    if (!ok)
    {
        throw std::runtime_error(message);
    }
}

void selftest()
{
    Gravity central(0);
    Vec r{3000, 4000, 5000};
    const double d = norm(r);
    const Vec expected = (-central.gm / (d * d * d)) * r;
    require(norm(central.acceleration(r) - expected) < 1e-17, "Central gravity test failed");
    Gravity g(2);
    g.set(2, 0, -2.2491270601602e-5);
    g.set(2, 2, 1.249e-5, -3.4e-6);
    for (Vec p : std::vector<Vec>{{3000, 4000, 5000}, {0, 0, 3000}, {0, 0, -3000}, {3000, 0, 0}})
    {
        Vec numerical{};
        for (int i = 0; i < 3; ++i)
        {
            Vec hi = p, lo = p;
            hi[i] += .001;
            lo[i] -= .001;
            numerical[i] = (g.potential(hi).v - g.potential(lo).v) / .002;
        }
        require(norm(g.acceleration(p) - numerical) < 1e-10, "Gravity derivative/pole test failed");
    }
    const auto oscillator = [](double, const State& s)
    { return State{s[3], s[4], s[5], -s[0], -s[1], -s[2]}; };
    Integrator<decltype(oscillator)> integ(oscillator, 0, {1, 0, 0, 0, 1, 0}, .1);
    const auto end = integ.to(2 * 3.14159265358979323846);
    require(std::abs(end[0] - 1) < 1e-8 && std::abs(end[1]) < 1e-8,
            "Integrator periodic orbit failed");
    std::cout << "PASS: central gravity, harmonic derivatives, both poles, RK8 integrator\n";
#if MESSENGER_STAGE >= 2
    CubicSpline line({0, 1, 2, 3}, {1, 3, 5, 7});
    require(std::abs(line(1.25) - 3.5) < 1e-14, "Spline linear reproduction failed");
#endif
}
#if MESSENGER_STAGE >= 2
void runLightTime(const std::string& data, const std::string& out, double epoch, double hours,
                  double interval, int station, const std::string& source, double integrationStep)
{
    Stations stations(data);
    SpacecraftSpline trajectory(epoch, epoch + hours * 3600);
    Gravity gravity(160);
    gravity.load(data + "/HgM009.sha");
    Dynamics dynamics(gravity);
    std::optional<NumericTrajectory> numeric;
    const double begin = epoch - 3600;
    if (source == "integrated")
    {
        numeric.emplace(dynamics, begin, epoch + hours * 3600, ephemeris(-236, begin),
                        integrationStep);
    }
    LightTime solver(stations, [&](Epoch et) { return numeric ? (*numeric)(et) : trajectory(et); });
    double orthogonality = 0;
    const auto site = stations.geocentric(station, epoch, &orthogonality);
    require(orthogonality < 1e-14 && norm(site) > 6300 && norm(site) < 6400,
            "Invalid Earth orientation/station units");

    std::filesystem::create_directories(out);
    std::ofstream f(out + "/lighttime.csv");
    f << std::setprecision(16)
      << "elapsed_s,up_s,down_s,roundtrip_s,up_shapiro_us,down_shapiro_us,up_closure_ns,down_"
         "closure_ns,up_iterations,down_iterations,reference_position_difference_m\n";
    double closure = 0, splineError = 0;
    for (double t = 0;; t = std::min(t + interval, hours * 3600))
    {
        const auto path = solver.solve(epoch + t, station, station);
        const double err =
            1000 * norm((numeric ? numeric->mercuryRelative(epoch + t - path.down_s)
                                 : trajectory.mercuryRelative(epoch + t - path.down_s)) -
                        pos(ephemeris(-236, epoch + t - path.down_s)));
        splineError = std::max(splineError, err);
        closure = std::max({closure, std::abs(path.up_closure_s), std::abs(path.down_closure_s)});
        f << t << ',' << path.up_s << ',' << path.down_s << ',' << path.up_s + path.down_s << ','
          << path.up_shapiro_s * 1e6 << ',' << path.down_shapiro_s * 1e6 << ','
          << path.up_closure_s * 1e9 << ',' << path.down_closure_s * 1e9 << ','
          << path.up_iterations << ',' << path.down_iterations << ',' << err << '\n';
        if (t >= hours * 3600)
        {
            break;
        }
    }
    require(closure < 1e-8, "Light-time closure exceeds tolerance");
    std::ofstream summary(out + "/summary.json");
    summary << std::setprecision(15) << "{\n  \"station\": " << station
            << ",\n  \"spacecraft_source\": \"" << source
            << "\",\n  \"eop_interpolation\": \"natural cubic spline\",\n  \"eop_units\": "
               "\"xp,yp,dX,dY: mas; TAI-UT1: ms\",\n  \"max_closure_ns\": "
            << closure * 1e9 << ",\n  \"max_reference_position_difference_m\": " << splineError
            << ",\n  \"rotation_orthogonality\": " << orthogonality << "\n}\n";
    std::cout << "Light-time complete: closure " << closure * 1e9 << " ns, reference difference "
              << splineError << " m\n";
}
#endif

// Настройки расчёта: длительность, шаг, степень гармоник и пути к файлам.
struct RunOptions
{
    std::string command = "orbit";
    std::string data = "data";
    std::string out;
    std::string trajectorySource = "integrated";
    double hours = 20;
    double interval = 60;
    double step = 60;
    double fitHours = 0;
    double countScale = .001;
    int degree = 160;
    int station = 34;
    int stride = 1;
    int points = 180;
    int iterations = 12;
};

// Чтение параметров запуска; без аргументов используются значения по умолчанию.
RunOptions parseArguments(int argc, char** argv)
{
    RunOptions options;
    if (argc > 1)
    {
        options.command = argv[1];
    }

    for (int i = 2; i < argc; ++i)
    {
        const std::string key = argv[i];
        if (i + 1 >= argc)
        {
            throw std::runtime_error("Missing option value");
        }
        const std::string value = argv[++i];

        if (key == "--data")
        {
            options.data = value;
        }
        else if (key == "--out")
        {
            options.out = value;
        }
        else if (key == "--hours")
        {
            options.hours = std::stod(value);
        }
        else if (key == "--interval")
        {
            options.interval = std::stod(value);
        }
        else if (key == "--degree")
        {
            options.degree = std::stoi(value);
        }
        else if (key == "--step")
        {
            options.step = std::stod(value);
        }
        else if (key == "--trajectory")
        {
            options.trajectorySource = value;
        }
        else if (key == "--count-scale")
        {
            options.countScale = std::stod(value);
        }
        else if (key == "--station")
        {
            options.station = std::stoi(value);
        }
        else if (key == "--stride")
        {
            options.stride = std::stoi(value);
        }
        else if (key == "--points")
        {
            options.points = std::stoi(value);
        }
        else if (key == "--iterations")
        {
            options.iterations = std::stoi(value);
        }
        else if (key == "--fit-hours")
        {
            options.fitHours = std::stod(value);
        }
        else
        {
            throw std::runtime_error("Unknown option: " + key);
        }
    }

    if (options.out.empty())
    {
        const std::string results = "results/corrected/";
        if (options.command == "fit")
        {
            options.out = results + "stage_1_4";
        }
        else if (options.command == "doppler")
        {
            options.out =
                results + (options.trajectorySource == "reference" ? "stage_1_3_reference"
                                                                   : "stage_1_3_integrated");
        }
        else if (options.command == "lighttime")
        {
            options.out = results + "stage_1_2";
        }
        else
        {
            options.out = results + "stage_1";
        }
    }

    if (options.command == "selftest")
    {
        return options;
    }

    require(options.command == "orbit" || options.command == "lighttime" ||
                options.command == "doppler" || options.command == "fit",
            "Available commands: orbit, lighttime, doppler, fit, selftest");

    require(std::isfinite(options.hours) && std::isfinite(options.interval) &&
                std::isfinite(options.step) && std::isfinite(options.fitHours) &&
                std::isfinite(options.countScale) && options.countScale > 0 && options.hours > 0 &&
                options.interval > 0 && options.step > 0 && options.degree >= 0 &&
                options.degree <= 160 && options.stride > 0 && options.points >= 12 &&
                options.iterations >= 1 && options.fitHours >= 0,
            "Invalid integration options");

    require(options.trajectorySource == "integrated" || options.trajectorySource == "reference",
            "Unknown trajectory source");
    require(options.command == "orbit" ||
                (MESSENGER_STAGE >= 2 && options.command == "lighttime") ||
                (MESSENGER_STAGE >= 3 && options.command == "doppler") ||
                (MESSENGER_STAGE >= 4 && options.command == "fit"),
            "Command unavailable in this cumulative version");
    require(options.command != "fit" || options.trajectorySource == "integrated",
            "Fitting requires the integrated trajectory");

    return options;
}

void runOrbit(const RunOptions& options)
{
    const std::string& data = options.data;
    const std::string& out = options.out;
    const double hours = options.hours;
    const double interval = options.interval;
    const double step = options.step;
    const int degree = options.degree;

    require(interval >= step && std::abs(interval / step - std::round(interval / step)) < 1e-10,
            "Output interval must be a multiple of the RK8 step");

    // Загрузка коэффициентов неоднородного поля Меркурия и создание модели движения.
    Gravity gravity(degree);
    gravity.load(data + "/HgM009.sha");
    Dynamics dynamics(gravity);

    // Начальный момент и шесть компонент состояния спутника из эфемерид.
    const double epoch = jdTDBtoET(2455646.83);
    const auto initial = ephemeris(-236, epoch);

    std::filesystem::create_directories(out);
    std::ofstream csv(out + "/orbit.csv");
    if (!csv)
    {
        throw std::runtime_error("Cannot write orbit.csv");
    }
    csv << "elapsed_s,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s,reference_x_km,reference_y_km,"
           "reference_z_km,position_error_m,velocity_error_m_s\n"
        << std::setprecision(16);

    // Интегратор хранит текущее состояние и продолжает расчёт с достигнутого времени.
    Integrator<Dynamics> integrator(dynamics, epoch, initial, step);
    double sumSquaredPositionErrorsM2 = 0, maxPositionErrorM = 0, positionErrorM = 0;
    size_t sampleCount = 0;

    auto startTime = std::chrono::steady_clock::now();
    for (double elapsedSeconds = 0;;
         elapsedSeconds = std::min(elapsedSeconds + interval, hours * 3600))
    {
        // Расчётное состояние и независимое эталонное состояние для того же момента.
        const auto state = integrator.to(epoch + elapsedSeconds),
                   reference = ephemeris(-236, epoch + elapsedSeconds);
        // Длина вектора отклонения в метрах, максимум и сумма квадратов для итогового СКО.
        positionErrorM = 1000 * norm(pos(state) - pos(reference));
        maxPositionErrorM = std::max(maxPositionErrorM, positionErrorM);
        sumSquaredPositionErrorsM2 += positionErrorM * positionErrorM;
        ++sampleCount;
        csv << elapsedSeconds;
        for (double x : state)
        {
            csv << ',' << x;
        }
        for (double x : pos(reference))
        {
            csv << ',' << x;
        }
        csv << ',' << positionErrorM << ',' << 1000 * norm(vel(state) - vel(reference)) << '\n';
        if (elapsedSeconds >= hours * 3600)
        {
            break;
        }
    }

    double runtimeSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
    // Сохранение настроек и итоговых показателей точности в отдельный файл.
    std::ofstream summary(out + "/summary.json");
    summary << std::setprecision(15)
            << "{\n  \"epoch_jd_tdb\": 2455646.83,\n  \"duration_hours\": " << hours
            << ",\n  \"gravity_degree\": " << degree
            << ",\n  \"solar_radiation_pressure\": false,\n  \"relativistic_dynamics\": \"EIH "
               "beta=gamma=1\",\n  \"integrator\": \"Dormand-Prince RK8 (13 stages)\",\n  "
               "\"step_seconds\": "
            << step << ",\n  \"GM_mercury_km3_s2\": " << gravity.gm
            << ",\n  \"mercury_orientation\": \"Genova 2019 pole, Stark 2015 rotation formula\",\n "
               " \"position_max_m\": "
            << maxPositionErrorM
            << ",\n  \"position_rms_m\": " << std::sqrt(sumSquaredPositionErrorsM2 / sampleCount)
            << ",\n  \"position_final_m\": " << positionErrorM
            << ",\n  \"accepted_steps\": " << integrator.stats.accepted
            << ",\n  \"rejected_steps\": " << integrator.stats.rejected
            << ",\n  \"runtime_s\": " << runtimeSeconds << "\n}\n";
    std::cout << std::setprecision(10) << "Degree " << degree << ", duration " << hours
              << " h, max error " << maxPositionErrorM << " m, RMS "
              << std::sqrt(sumSquaredPositionErrorsM2 / sampleCount) << " m, final "
              << positionErrorM << " m\n"
              << "Steps: " << integrator.stats.accepted << " accepted, "
              << integrator.stats.rejected << " rejected; runtime " << runtimeSeconds << " s\n";
}

int main(int argc, char** argv)
{
    try
    {
        const RunOptions options = parseArguments(argc, argv);

        if (options.command == "selftest")
        {
            selftest();
            return 0;
        }

        // Загрузка эфемерид, расчёт орбиты и освобождение ресурсов библиотеки.
        loadKernels(options.data);

#if MESSENGER_STAGE >= 2
        if (options.command == "lighttime")
        {
            runLightTime(options.data, options.out, jdTDBtoET(2455646.83), options.hours,
                         options.interval, options.station, options.trajectorySource, options.step);
            kclear_c();
            return 0;
        }
#endif

#if MESSENGER_STAGE >= 3
        if (options.command == "doppler")
        {
            runDoppler(options.data, options.out, options.stride, options.countScale,
                       options.trajectorySource, options.step);
            kclear_c();
            return 0;
        }
#endif

#if MESSENGER_STAGE >= 4
        if (options.command == "fit")
        {
            runFit(options.data, options.out, options.points, options.iterations, options.fitHours,
                   options.countScale, options.step);
            kclear_c();
            return 0;
        }
#endif

        runOrbit(options);
        kclear_c();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}

