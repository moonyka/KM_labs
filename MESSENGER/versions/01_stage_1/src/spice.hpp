#pragma once
#include "gravity.hpp"
#include "eih.hpp"
#include "epoch.hpp"
#include "mercury.hpp"
#include "shgrav.hpp"
#include <SpiceUsr.h>
#include <filesystem>
#include <string>
#include <map>
#include <fstream>

inline void spiceCheck()
{
    if (failed_c())
    {
        char msg[1841];
        getmsg_c("LONG", 1841, msg);
        reset_c();
        throw std::runtime_error(msg);
    }
}

inline void loadKernel(const std::string& path)
{
    if (!std::filesystem::exists(path))
    {
        throw std::runtime_error("Missing kernel: " + path);
    }
    furnsh_c(path.c_str());
    spiceCheck();
}

// Подключение файлов времени, параметров тел и эфемерид.
inline void loadKernels(const std::string& data)
{
    erract_c("SET", 0, const_cast<char*>("RETURN"));
    errprt_c("SET", 0, const_cast<char*>("NONE"));
    for (const auto& f : {"naif0012.tls", "pck00010.tpc", "de442t.bsp", "msgr_final.bsp"})
    {
        loadKernel(data + "/" + f);
    }
}

// Геометрическое состояние в осях J2000; по умолчанию начало координат — Меркурий.
inline State ephemeris(int target, double et, int observer = 199)
{
    State s{};
    double lt;
    spkez_c(target, et, "J2000", "NONE", observer, s.data(), &lt);
    spiceCheck();
    return s;
}

inline State ephemeris(int target, const Epoch& et, int observer = 199)
{
    State s = ephemeris(target, et.high, observer);
    for (int j = 0; j < 3; ++j)
    {
        s[j] += et.low * s[j + 3];
    }
    return s;
}

struct Body
{
    int id;
    double gm;
};

inline const std::vector<Body>& massiveBodies()
{
    static const std::vector<Body> values{
        {10, 132712440041.279419}, {199, 22031.868551}, {2, 324858.592000},    {399, 398600.435507},
        {301, 4902.800118},        {4, 42828.375816},   {5, 126712764.100000}, {6, 37940584.841800},
        {7, 5794556.400000},       {8, 6836527.100580}, {9, 975.500000}};
    return values;
}

inline double bodyGM(int id)
{
    for (const auto& b : massiveBodies())
    {
        if (b.id == id)
        {
            return b.gm;
        }
    }
    throw std::runtime_error("Unknown gravitational parameter");
}

inline double jdTDBtoET(double jd)
{
    return (jd - 2451545.) * 86400.;
}

// Из PCK берём измеренные параметры полюса, вращения и основной либрации.
// Готовые высшие амплитуды и отдельные фазы из PCK здесь не используются.
inline MercuryParameters loadMercuryParameters(const std::filesystem::path& orbitInput)
{
    MercuryParameters model;
    const auto readPool = [](const char* name, int count, double* values)
    {
        SpiceInt received = 0;
        SpiceBoolean found = SPICEFALSE;
        gdpool_c(name, 0, count, &received, values, &found);
        spiceCheck();
        if (!found || received != count)
        {
            throw std::runtime_error(std::string("Missing PCK parameter: ") + name);
        }
        for (int i = 0; i < count; ++i)
        {
            if (!std::isfinite(values[i]))
            {
                throw std::runtime_error(std::string("Invalid PCK parameter: ") + name);
            }
        }
    };
    readPool("BODY199_POLE_RA", 3, model.poleRA.data());
    readPool("BODY199_POLE_DEC", 3, model.poleDec.data());
    readPool("BODY199_PM", 3, model.primeMeridian.data());
    readPool("BODY199_NUT_PREC_PM", 1, &model.fundamentalLibrationDeg);
    double argument[2];
    readPool("BODY1_NUT_PREC_ANGLES", 2, argument);
    model.phaseAtJ2000Deg = argument[0];
    model.meanMotionDegDay = argument[1] / 36525.;

    // Эксцентриситет — отдельное исходное наблюдаемое значение из статьи Stark.
    std::ifstream input(orbitInput);
    if (!input)
    {
        throw std::runtime_error("Cannot open Mercury orbit inputs: " + orbitInput.string());
    }
    bool hasEccentricity = false;
    std::string key;
    while (input >> key)
    {
        if (key.front() == '#')
        {
            std::string comment;
            std::getline(input, comment);
            continue;
        }
        if (key != "eccentricity" || hasEccentricity || !(input >> model.eccentricity))
        {
            throw std::runtime_error("Invalid Mercury orbit inputs");
        }
        hasEccentricity = true;
    }
    if (!hasEccentricity)
    {
        throw std::runtime_error("Missing Mercury eccentricity");
    }
    // Stark (3): все пять отношений вычисляются из e один раз при создании модели.
    model.relativeAmplitudes = calculateLibrationRatios(model.eccentricity);
    return model;
}

class Dynamics
{
    const Gravity& gravity_;
    ShGrav harmonics_;
    MercuryParameters rotationParameters_;

    struct Environment
    {
        State mercury;
        std::vector<GravBody> bodies;
        Vec mercuryCorrection;
    };

    mutable std::map<double, Environment> cache_;

    // Положения и скорости массивных тел относительно барицентра для расчёта возмущений.
    const Environment& environmentAt(double et) const
    {
        const auto found = cache_.find(et);
        if (found != cache_.end())
        {
            return found->second;
        }
        Environment e{};
        e.mercury = ephemeris(199, et, 0);
        for (const auto& b : massiveBodies())
        {
            const auto s = b.id == 199 ? e.mercury : ephemeris(b.id, et, 0);
            e.bodies.push_back({b.id, b.gm, pos(s), vel(s)});
        }
        // Мойер (4-26), (4-27): подготовка U_j и ньютоновских a_j для этой эпохи.
        // Повторный запрос той же эпохи использует кеш; состояние КА не кешируется.
        prepareEIH(e.bodies);
        e.mercuryCorrection = eihCorrection(pos(e.mercury), vel(e.mercury), e.bodies, 199);
        if (cache_.size() > 50000)
        {
            cache_.clear();
        }
        return cache_.emplace(et, std::move(e)).first->second;
    }

public:
    explicit Dynamics(const Gravity& g)
        : gravity_(g),
          harmonics_(shGravCreate(g.degree(),
                                  g.gm * DAY_SECONDS * DAY_SECONDS / (AU_KM * AU_KM * AU_KM),
                                  g.radius / AU_KM)),
          rotationParameters_(loadMercuryParameters(
              std::filesystem::path(g.coefficientFile).parent_path().parent_path() /
              "config" / "mercury_orbit.txt"))
    {
        if (g.coefficientFile.empty())
        {
            throw std::runtime_error("Dynamics requires loaded harmonic coefficients");
        }
        if (shGravLoadCoefficients(harmonics_, g.coefficientFile.c_str()) < 0)
        {
            throw std::runtime_error("Cannot load harmonic coefficients");
        }
    }

    Vec acceleration(double et, const State& s) const
    {
        // Поле Меркурия: центральное притяжение и гармоники.
        const Vec r = pos(s);
        double rotation[3][3];
        mercuryRotation(et, rotation, rotationParameters_);
        const double distance = norm(r);
        // Поворот во вращающиеся оси Меркурия и перевод километров в астрономические единицы.
        const Vec bodyFixed = (1 / AU_KM) * mv(rotation, r);
        double harmonicAcceleration[3];
        shGravCalculate(harmonics_, bodyFixed[0], bodyFixed[1], bodyFixed[2], harmonicAcceleration);
        // Перевод ускорения в км/с²; ниже оно поворачивается обратно в оси J2000.
        const Vec harmonicInKm =
            (AU_KM / (DAY_SECONDS * DAY_SECONDS)) *
            Vec{harmonicAcceleration[0], harmonicAcceleration[1], harmonicAcceleration[2]};
        Vec a = (-gravity_.gm / (distance * distance * distance)) * r + mtv(rotation, harmonicInKm);
        // Дифференциальное ньютоновское притяжение внешних тел.
        const auto& e = environmentAt(et);
        const Vec mercury = pos(e.mercury);
        for (const auto& b : e.bodies)
        {
            if (b.id != 199)
            {
                const Vec rb = b.r - mercury, d = rb - r;
                const double rn = norm(rb), dn = norm(d);
                a = a + b.gm * ((1 / (dn * dn * dn)) * d - (1 / (rn * rn * rn)) * rb);
            }
        }
        // Разность релятивистских поправок аппарата и Меркурия.
        return a + eihCorrection(mercury + r, vel(e.mercury) + vel(s), e.bodies) -
               e.mercuryCorrection;
    }

    // Правая часть: производные координат равны скоростям, производные скоростей — ускорениям.
    State operator()(double time, const State& state) const
    {
        const Vec totalAcceleration = acceleration(time, state);
        return {
            state[3],
            state[4],
            state[5], // Производная координат — скорость.
            totalAcceleration[0],
            totalAcceleration[1],
            totalAcceleration[2] // Производная скорости — ускорение.
        };
    }
};

