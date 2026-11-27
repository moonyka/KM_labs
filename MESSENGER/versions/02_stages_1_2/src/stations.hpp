#pragma once
#include "spice.hpp"
#include "spline.hpp"
#include <sofa.h>
#include <map>
#include <sstream>
#include <fstream>

class Stations
{
    struct Station
    {
        Vec r, v;
    };

    std::map<int, Station> stations_;
    std::array<CubicSpline, 5> eop_;

public:
    explicit Stations(const std::string& data)
    {
        std::ifstream f(data + "/stations_19jul2013.txt");
        if (!f)
        {
            throw std::runtime_error("Missing station coordinates");
        }
        std::string line;
        while (std::getline(f, line))
        {
            std::istringstream row(line);
            std::string id, diam, type;
            Station s{};
            if (!(row >> id >> diam >> type >> s.r[0] >> s.r[1] >> s.r[2] >> s.v[0] >> s.v[1] >>
                  s.v[2]))
            {
                continue;
            }
            if (id.find('*') != std::string::npos)
            {
                continue;
            }
            try
            {
                stations_[std::stoi(id)] = s;
            }
            catch (const std::invalid_argument&)
            {
            }
        }
        std::ifstream e(data + "/eop2.txt");
        size_t n;
        e >> n;
        if (!e || n < 3)
        {
            throw std::runtime_error("Invalid EOP table");
        }
        std::vector<double> x(n);
        std::array<std::vector<double>, 5> y;
        for (auto& a : y)
        {
            a.resize(n);
        }
        for (size_t i = 0; i < n; ++i)
        {
            e >> x[i];
            for (auto& a : y)
            {
                e >> a[i];
            }
            if (!e)
            {
                throw std::runtime_error("Truncated EOP table");
            }
        }
        for (int j = 0; j < 5; ++j)
        {
            eop_[j] = CubicSpline(x, y[j]);
        }
    }

    Vec geocentric(int id, double et, double* orthogonality = nullptr) const
    {
        const auto it = stations_.find(id);
        if (it == stations_.end())
        {
            throw std::runtime_error("Unknown DSN station " + std::to_string(id));
        }
        // Письмо от 28.09: TDB -> TT из de442t, затем TT -> TAI через SOFA.
        const double ttMinusTdb = ephemeris(1000000001, et, 1000000000)[0];
        const double tt1 = 2451545., tt2 = (et + ttMinusTdb) / 86400.;
        double tai1, tai2;
        if (iauTttai(tt1, tt2, &tai1, &tai2) != 0)
        {
            throw std::runtime_error("SOFA time conversion failed");
        }
        const double mjd = (tai1 - 2400000.5) + tai2;
        const double mas2rad = 3.14159265358979323846 / (180 * 3600 * 1000.);
        const double xp = eop_[0](mjd) * mas2rad, yp = eop_[1](mjd) * mas2rad;
        const double taiMinusUt1 = eop_[2](mjd) * .001;
        const double dx = eop_[3](mjd) * mas2rad, dy = eop_[4](mjd) * mas2rad;

        // CIP по iauXy06; поправки dX/dY из файла преподавателя заданы в mas.
        double x, y;
        iauXy06(tt1, tt2, &x, &y);
        x += dx;
        y += dy;
        const double ut11 = tai1, ut12 = tai2 - taiMinusUt1 / 86400.;
        double rc2t[3][3];
        // C2T из письма; ниже транспонирование даёт нужное направление T2C.
        iauC2txy(tt1, tt2, ut11, ut12, x, y, xp, yp, rc2t);
        if (orthogonality)
        {
            *orthogonality = 0;
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    double v = 0;
                    for (int k = 0; k < 3; ++k)
                    {
                        v += rc2t[k][i] * rc2t[k][j];
                    }
                    *orthogonality = std::max(*orthogonality, std::abs(v - (i == j ? 1 : 0)));
                }
            }
        }
        // Файл станций: метры и метры/год, поэтому интервал переводим в годы.
        const double years = ((tt1 - 2453371.5) + tt2) / 365.25;
        Vec itrf = (.001) * (it->second.r + years * it->second.v);
        return mtv(rc2t, itrf);
    }

    Vec geocentric(int id, const Epoch& et) const
    {
        return geocentric(id, et.high) +
               (.5 * et.low) * (geocentric(id, et.high + 1) - geocentric(id, et.high - 1));
    }

    Vec barycentric(int id, const Epoch& et) const
    {
        return pos(ephemeris(399, et, 0)) + geocentric(id, et);
    }

    Vec barycentric(int id, double et) const
    {
        return barycentric(id, Epoch(et));
    }
};
