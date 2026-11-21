#include "../src/trajectory.hpp"
#include "../src/lighttime.hpp"
#include <iostream>
void verify(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
int main(int argc,char**argv){try{
    const std::string data=argc>1?argv[1]:"data";loadKernels(data);
    double ra[3]={281.0082,-.032808,0},dec[3]={61.4164,-.0048464,0},pm[3]={329.5469,6.1385025,0};
    double phase[10],libration[5];const double amplitude[5]={1.,-.105231253725171,-.0103439039634868,-.0023786959435392,-.000535306724141346};
    for(int k=1;k<=5;++k){phase[2*k-2]=std::fmod(k*4.09233445*42.71182,360.);phase[2*k-1]=k*4.09233445*36525.;libration[k-1]=(40./3600.)*amplitude[k-1];}
    pdpool_c("BODY199_POLE_RA",3,ra);pdpool_c("BODY199_POLE_DEC",3,dec);pdpool_c("BODY199_PM",3,pm);pdpool_c("BODY199_NUT_PREC_PM",5,libration);pdpool_c("BODY1_NUT_PREC_ANGLES",10,phase);spiceCheck();
    double rotationError=0;
    for(double t:{0.,354358512.,-1e8}){double own[3][3],spice[3][3];mercuryRotation(t,own);pxform_c("J2000","IAU_MERCURY",t,spice);spiceCheck();for(int i=0;i<3;++i)for(int j=0;j<3;++j)rotationError=std::max(rotationError,std::abs(own[i][j]-spice[i][j]));}
    verify(rotationError<2e-12,"Mercury Euler rotation differs from independent CSPICE PCK evaluation");
    Gravity gravity(160);gravity.load(data+"/HgM009.sha");Dynamics dynamics(gravity);const double begin=jdTDBtoET(2455646.83);const State initial=ephemeris(-236,begin);
    double pericenter=begin,minRadius=norm(pos(initial));for(int j=1;j<=1200;++j){const double t=begin+60*j,r=norm(pos(ephemeris(-236,t)));if(r<minRadius){minRadius=r;pericenter=t;}}
    const double arcBegin=pericenter-300;const State arcInitial=ephemeris(-236,arcBegin);
    NumericTrajectory dense(dynamics,arcBegin,arcBegin+600,arcInitial);Integrator<Dynamics> direct(dynamics,arcBegin,arcInitial,5);double interpolationError=0;
    for(int i=0;i<20;++i){const double t=arcBegin+15+30*i;interpolationError=std::max(interpolationError,1000*norm(dense.mercuryRelative(t)-pos(direct.to(t))));}
    verify(interpolationError<.02,"Quintic trajectory interpolation exceeds 2 cm on validation arc");
    Stations stations(data);SpacecraftSpline reference(begin,begin+600);LightTime light(stations,[&](Epoch t){return reference(t);});const auto path=light.solve(Epoch(begin,.00000001),34,34);
    verify(std::max(std::abs(path.up_closure_s),std::abs(path.down_closure_s))<1e-8,"Light-time equation does not close");
    std::cout<<"PASS: independent PCK rotation error "<<rotationError<<", dense interpolation "<<interpolationError<<" m, split-epoch two-way light time\n";kclear_c();return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
