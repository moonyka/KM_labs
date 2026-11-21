#pragma once
#include "stations.hpp"
#include <functional>
inline constexpr double LIGHT_SPEED=299792.458;
// Natural cubic splines of the Mercury-relative reference trajectory. Padding
// keeps natural end conditions away from the requested observation interval.
class SpacecraftSpline {
    double base_;std::array<CubicSpline,3> coordinates_;
public:
    SpacecraftSpline(double from,double to,double step=10):base_(from-3600){
        std::vector<double> t;std::array<std::vector<double>,3> r;
        for(double et=base_;et<=to+120+step;et+=step){const auto state=ephemeris(-236,et);t.push_back(et-base_);for(int j=0;j<3;++j)r[j].push_back(state[j]);}
        for(int j=0;j<3;++j)coordinates_[j]=CubicSpline(t,r[j]);
    }
    Vec mercuryRelative(const Epoch&et)const{const double t=(et.high-base_)+et.low;return {coordinates_[0](t),coordinates_[1](t),coordinates_[2](t)};}
    Vec mercuryRelative(double et)const{return mercuryRelative(Epoch(et));}
    Vec operator()(const Epoch&et)const{return pos(ephemeris(199,et,0))+mercuryRelative(et);}
};
struct LightPath {Epoch receive;double receive_et,up_s,down_s,up_shapiro_s,down_shapiro_s,up_closure_s,down_closure_s;int up_iterations,down_iterations;Epoch transmitTime()const{return receive.shifted(-down_s-up_s);}Epoch bounceTime()const{return receive.shifted(-down_s);}double transmitET()const{return transmitTime().high;}double bounceET()const{return bounceTime().high;}};
class LightTime {
    const Stations& stations_;std::function<Vec(Epoch)> spacecraft_;std::vector<Body> bodies_;
public:
    LightTime(const Stations&s,std::function<Vec(Epoch)>sc):stations_(s),spacecraft_(std::move(sc)){for(int id:{10,199,2,399,301,4,5,6,7,8,9})bodies_.push_back({id,bodyGM(id)});}
    // Moyer (2000), Eq. 8-55, gamma=1. Solar term includes bending correction.
    double shapiro(const Vec&r1,Epoch t1,const Vec&r2,Epoch t2)const{
        double delay=0;
        for(const auto&b:bodies_){const Vec p=r1-pos(ephemeris(b.id,t1,0)),q=r2-pos(ephemeris(b.id,t2,0));const double sum=norm(p)+norm(q),length=norm(q-p),correction=b.id==10?2*b.gm/(LIGHT_SPEED*LIGHT_SPEED):0;
            if(sum-length+correction<=0)throw std::runtime_error("Invalid relativistic light path");delay+=2*b.gm/std::pow(LIGHT_SPEED,3)*std::log((sum+length+correction)/(sum-length+correction));}
        return delay;
    }
    LightPath solve(Epoch t3,int transmitStation,int receiveStation)const{
        LightPath p{};p.receive=t3;p.receive_et=t3.high;const Vec rx=stations_.barycentric(receiveStation,t3);double tau=norm(rx-spacecraft_(t3))/LIGHT_SPEED;
        for(int i=0;i<20;++i){const Epoch t2=t3.shifted(-tau);const Vec sc=spacecraft_(t2);const double sr=shapiro(sc,t2,rx,t3),next=norm(rx-sc)/LIGHT_SPEED+sr;p.down_iterations=i+1;if(std::abs(next-tau)<1e-10){tau=next;break;}tau=next;if(i==19)throw std::runtime_error("Downlink did not converge");}
        p.down_s=tau;const Epoch t2=t3.shifted(-tau);const Vec sc=spacecraft_(t2);p.down_shapiro_s=shapiro(sc,t2,rx,t3);p.down_closure_s=tau-norm(rx-sc)/LIGHT_SPEED-p.down_shapiro_s;
        tau=norm(sc-stations_.barycentric(transmitStation,t2))/LIGHT_SPEED;
        for(int i=0;i<20;++i){const Epoch t1=t2.shifted(-tau);const Vec tx=stations_.barycentric(transmitStation,t1);const double sr=shapiro(tx,t1,sc,t2),next=norm(sc-tx)/LIGHT_SPEED+sr;p.up_iterations=i+1;if(std::abs(next-tau)<1e-10){tau=next;break;}tau=next;if(i==19)throw std::runtime_error("Uplink did not converge");}
        p.up_s=tau;const Epoch t1=t2.shifted(-tau);const Vec tx=stations_.barycentric(transmitStation,t1);p.up_shapiro_s=shapiro(tx,t1,sc,t2);p.up_closure_s=tau-norm(sc-tx)/LIGHT_SPEED-p.up_shapiro_s;return p;
    }
};
