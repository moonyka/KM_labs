#pragma once
#include "spice.hpp"
#include "spline.hpp"
#include <sofa.h>
#include <map>
#include <sstream>
#include <fstream>

class Stations {
    struct Station{Vec r,v;};std::map<int,Station> stations_;std::array<CubicSpline,5> eop_;
public:
    explicit Stations(const std::string&data){
        std::ifstream f(data+"/stations_19jul2013.txt");if(!f)throw std::runtime_error("Missing station coordinates");std::string line;
        while(std::getline(f,line)){std::istringstream row(line);std::string id,diam,type;Station s{};if(!(row>>id>>diam>>type>>s.r[0]>>s.r[1]>>s.r[2]>>s.v[0]>>s.v[1]>>s.v[2]))continue;if(id.find('*')!=std::string::npos)continue;try{stations_[std::stoi(id)]=s;}catch(const std::invalid_argument&){} }
        std::ifstream e(data+"/eop2.txt");size_t n;e>>n;if(!e||n<3)throw std::runtime_error("Invalid EOP table");std::vector<double>x(n);std::array<std::vector<double>,5>y;for(auto&a:y)a.resize(n);
        for(size_t i=0;i<n;++i){e>>x[i];for(auto&a:y)e>>a[i];if(!e)throw std::runtime_error("Truncated EOP table");}for(int j=0;j<5;++j)eop_[j]=CubicSpline(x,y[j]);
    }
    Vec geocentric(int id,double et,double*orthogonality=nullptr)const{
        const auto it=stations_.find(id);if(it==stations_.end())throw std::runtime_error("Unknown DSN station "+std::to_string(id));
        // de442t contains TT-TDB; EOP epochs in the supplied assignment are TAI.
        const double ttMinusTdb=ephemeris(1000000001,et,1000000000)[0];
        const double tai=et+ttMinusTdb-32.184;double tt1,tt2;
        if(iauTaitt(2451545.,tai/86400.,&tt1,&tt2)!=0)throw std::runtime_error("SOFA time conversion failed");
        const double mjd=(2451545.-2400000.5)+tai/86400.,mas2rad=3.14159265358979323846/(180*3600*1000.);
        const double xp=eop_[0](mjd)*mas2rad,yp=eop_[1](mjd)*mas2rad,tai_ut1=eop_[2](mjd)*.001,dx=eop_[3](mjd)*mas2rad,dy=eop_[4](mjd)*mas2rad;
        double x,y,s,rc2i[3][3],rpom[3][3],rc2t[3][3];iauXys06a(tt1,tt2,&x,&y,&s);iauC2ixys(x+dx,y+dy,s,rc2i);
        iauPom00(xp,yp,iauSp00(tt1,tt2),rpom);const double era=iauEra00(2451545.,(tai-tai_ut1)/86400.);iauC2tcio(rc2i,era,rpom,rc2t);
        if(orthogonality){*orthogonality=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j){double v=0;for(int k=0;k<3;++k)v+=rc2t[k][i]*rc2t[k][j];*orthogonality=std::max(*orthogonality,std::abs(v-(i==j?1:0)));}}
        const double years=((tt1-2453371.5)+tt2)/365.25;Vec itrf=(.001)*(it->second.r+years*it->second.v);return mtv(rc2t,itrf);
    }
    Vec barycentric(int id,double et)const{return pos(ephemeris(399,et,0))+geocentric(id,et);}
};
