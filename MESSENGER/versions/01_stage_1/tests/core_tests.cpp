#include "../src/gravity.hpp"
#include "../src/shgrav.hpp"
#include "../src/eih.hpp"
#include "../src/mercury.hpp"
#include "../src/integrator.hpp"
#include "../src/epoch.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
void check(bool condition,const char*message){if(!condition)throw std::runtime_error(message);}
int main(int argc,char**argv){try{
    const Epoch e(354067200.,.00000001);check(std::abs(difference(e,Epoch(354067200.))-1e-8)<1e-22,"Split epoch lost nanoseconds");
    check(std::abs(difference(e.shifted(.5),e)-.5)<1e-15,"Split epoch lost count duration");
    const auto exponential=[](double,const State&s){State d{};d[0]=s[0];return d;};
    const auto error=[&](double h){Integrator<decltype(exponential)> i(exponential,0,{1,0,0,0,0,0},h);return std::abs(i.to(5)[0]-std::exp(5.));};
    const double coarse=error(1),fine=error(.5);check(coarse>80*fine&&fine<1e-5,"RK8 order/convergence failed");
    const auto oscillator=[](double,const State&s){return State{s[3],s[4],s[5],-s[0],-s[1],-s[2]};};
    Integrator<decltype(oscillator)> orbit(oscillator,0,{1,0,0,0,1,0},.1);const auto end=orbit.to(2*3.14159265358979323846);
    check(std::abs(end[0]-1)<1e-9&&std::abs(end[1])<1e-9,"RK8 oscillator failed");check(std::abs(orbit.to(0)[0]-1)<1e-9,"RK8 backward integration failed");
    std::vector<GravBody> sun{{10,132712440041.279419,{0,0,0},{0,0,0}}};prepareEIH(sun);
    for(const auto&r:std::vector<Vec>{{5e7,0,0},{3e7,4e7,1e7}})for(const auto&v:std::vector<Vec>{{0,47,0},{10,-20,5}}){
        const double radius=norm(r),mu=sun[0].gm;const Vec expected=(mu/(LIGHT_SPEED_KM_S*LIGHT_SPEED_KM_S*radius*radius*radius))*((4*mu/radius-dot(v,v))*r+4*dot(r,v)*v);
        check(norm(eihCorrection(r,v,sun)-expected)<1e-24,"EIH Schwarzschild limit failed");
    }
    for(double t:{0.,354358512.,-1e8}){double R[3][3];mercuryRotation(t,R);for(int i=0;i<3;++i)for(int j=0;j<3;++j){double d=0;for(int k=0;k<3;++k)d+=R[i][k]*R[j][k];check(std::abs(d-(i==j?1.:0.))<1e-14,"Mercury rotation is not orthogonal");}
        const Vec north{R[2][0],R[2][1],R[2][2]};check(norm(mv(R,north)-Vec{0,0,1})<1e-14,"Mercury pole mapping failed");}
    if(argc>1){Gravity gravity(160);gravity.load(argv[1]);std::ostream*output=&std::cout;std::ofstream file;if(argc>2){file.open(argv[2]);output=&file;}
        *output<<std::setprecision(17)<<"x_km,y_km,z_km,ax_km_s2,ay_km_s2,az_km_s2\n";
        auto field=shGravCreate(160,gravity.gm*DAY_SECONDS*DAY_SECONDS/(AU_KM*AU_KM*AU_KM),2440./AU_KM);shGravLoadCoefficients(field,argv[1]);
        std::vector<Vec> probes{{3000,4000,5000},{2500,100,200},{-2700,150,-400},{0,0,3000},{0,0,-3000},{0.001,0,3000},{3000,0,0},{0,3000,0},{-2441,3,2}};
        for(int j=0;j<40;++j){const double latitude=(-85.+170.*j/39)*DEG_TO_RAD,longitude=137.5*j*DEG_TO_RAD,radius=2441.+400*j;probes.push_back({radius*std::cos(latitude)*std::cos(longitude),radius*std::cos(latitude)*std::sin(longitude),radius*std::sin(latitude)});}
        for(Vec r:probes){const Vec a=gravity.acceleration(r);check(norm(a-gravity.accelerationCartesian(r))<2e-16,"Degree-160 Legendre and Cartesian gradients disagree");
            double value[3];shGravCalculate(field,r[0]/AU_KM,r[1]/AU_KM,r[2]/AU_KM,value);const double distance=norm(r);const Vec api=(AU_KM/(DAY_SECONDS*DAY_SECONDS))*Vec{value[0],value[1],value[2]}+(-gravity.gm/(distance*distance*distance))*r;
            check(norm(api-a)<2e-16,"shGrav AU/day conversion or central-term accounting failed");for(double x:r)*output<<x<<',';*output<<a[0]<<','<<a[1]<<','<<a[2]<<'\n';}
    }
    std::cout<<"PASS: RK8 order, backward orbit, EIH Schwarzschild limit, Mercury rotation\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
