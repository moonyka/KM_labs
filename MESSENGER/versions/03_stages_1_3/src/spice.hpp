#pragma once
#include "gravity.hpp"
#include <SpiceUsr.h>
#include <filesystem>
#include <string>
#include <vector>
inline void spiceCheck(){if(failed_c()){char msg[1841];getmsg_c("LONG",1841,msg);reset_c();throw std::runtime_error(msg);}}
inline void loadKernel(const std::string& path){if(!std::filesystem::exists(path))throw std::runtime_error("Missing kernel: "+path);furnsh_c(path.c_str());spiceCheck();}
inline void loadKernels(const std::string& data){
    erract_c("SET",0,const_cast<char*>("RETURN"));errprt_c("SET",0,const_cast<char*>("NONE"));
    for(const auto&f:{"naif0012.tls","pck00010.tpc","gm_de431.tpc","de442t.bsp","msgr_final.bsp"})loadKernel(data+"/"+f);
}
inline State ephemeris(int target,double et,int observer=199){State s{};double lt;spkez_c(target,et,"J2000","NONE",observer,s.data(),&lt);spiceCheck();return s;}
inline double bodyGM(int id){double gm;SpiceInt n;bodvcd_c(id,"GM",1,&n,&gm);spiceCheck();return gm;}
inline double jdTDBtoET(double jd){return (jd-2451545.)*86400.;}
struct Body {int id;double gm;};
class Dynamics {
    const Gravity& gravity_;
    std::vector<Body> bodies_;
public:
    explicit Dynamics(const Gravity&g):gravity_(g){for(int id:{10,2,399,301,4,5,6,7,8})bodies_.push_back({id,bodyGM(id)});}
    Vec acceleration(double et,const State&s)const{
        const Vec r=pos(s);double rot[3][3];pxform_c("J2000","IAU_MERCURY",et,rot);spiceCheck();
        Vec a=mtv(rot,gravity_.acceleration(mv(rot,r)));
        for(const auto&b:bodies_){const Vec rb=pos(ephemeris(b.id,et)),d=rb-r;const double rn=norm(rb),dn=norm(d);a=a+b.gm*((1/(dn*dn*dn))*d-(1/(rn*rn*rn))*rb);}
        return a;
    }
    State operator()(double et,const State&s)const{State d{};const auto a=acceleration(et,s);for(int i=0;i<3;++i){d[i]=s[i+3];d[i+3]=a[i];}return d;}
};
