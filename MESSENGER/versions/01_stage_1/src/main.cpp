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

void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void selftest(){
    Gravity central(0);Vec r{3000,4000,5000};const double d=norm(r);const Vec expected=(-central.gm/(d*d*d))*r;
    require(norm(central.acceleration(r)-expected)<1e-17,"Central gravity test failed");
    Gravity g(2);g.set(2,0,-2.2491270601602e-5);g.set(2,2,1.249e-5,-3.4e-6);
    for(Vec p:std::vector<Vec>{{3000,4000,5000},{0,0,3000},{0,0,-3000},{3000,0,0}}){Vec numerical{};for(int i=0;i<3;++i){Vec hi=p,lo=p;hi[i]+=.001;lo[i]-=.001;numerical[i]=(g.potential(hi).v-g.potential(lo).v)/.002;}require(norm(g.acceleration(p)-numerical)<1e-10,"Gravity derivative/pole test failed");}
    const auto oscillator=[](double,const State&s){return State{s[3],s[4],s[5],-s[0],-s[1],-s[2]};};
    Integrator<decltype(oscillator)> integ(oscillator,0,{1,0,0,0,1,0},.1);const auto end=integ.to(2*3.14159265358979323846);
    require(std::abs(end[0]-1)<1e-8&&std::abs(end[1])<1e-8,"Integrator periodic orbit failed");
    std::cout<<"PASS: central gravity, harmonic derivatives, both poles, RK8 integrator\n";
#if MESSENGER_STAGE >= 2
    CubicSpline line({0,1,2,3},{1,3,5,7});require(std::abs(line(1.25)-3.5)<1e-14,"Spline linear reproduction failed");
#endif
}
#if MESSENGER_STAGE >= 2
void runLightTime(const std::string&data,const std::string&out,double epoch,double hours,double interval,int station,const std::string&source,double integrationStep){
    Stations stations(data);SpacecraftSpline trajectory(epoch,epoch+hours*3600);Gravity gravity(160);gravity.load(data+"/HgM009.sha");Dynamics dynamics(gravity);std::optional<NumericTrajectory> numeric;const double begin=epoch-3600;if(source=="integrated")numeric.emplace(dynamics,begin,epoch+hours*3600,ephemeris(-236,begin),integrationStep);LightTime solver(stations,[&](Epoch et){return numeric?(*numeric)(et):trajectory(et);});
    double orthogonality=0;const auto site=stations.geocentric(station,epoch,&orthogonality);require(orthogonality<1e-14&&norm(site)>6300&&norm(site)<6400,"Invalid Earth orientation/station units");
    std::filesystem::create_directories(out);std::ofstream f(out+"/lighttime.csv");f<<std::setprecision(16)<<"elapsed_s,up_s,down_s,roundtrip_s,up_shapiro_us,down_shapiro_us,up_closure_ns,down_closure_ns,up_iterations,down_iterations,reference_position_difference_m\n";double closure=0,splineError=0;
    for(double t=0;;t=std::min(t+interval,hours*3600)){const auto path=solver.solve(epoch+t,station,station);const double err=1000*norm((numeric?numeric->mercuryRelative(epoch+t-path.down_s):trajectory.mercuryRelative(epoch+t-path.down_s))-pos(ephemeris(-236,epoch+t-path.down_s)));splineError=std::max(splineError,err);closure=std::max({closure,std::abs(path.up_closure_s),std::abs(path.down_closure_s)});f<<t<<','<<path.up_s<<','<<path.down_s<<','<<path.up_s+path.down_s<<','<<path.up_shapiro_s*1e6<<','<<path.down_shapiro_s*1e6<<','<<path.up_closure_s*1e9<<','<<path.down_closure_s*1e9<<','<<path.up_iterations<<','<<path.down_iterations<<','<<err<<'\n';if(t>=hours*3600)break;}
    require(closure<1e-8,"Light-time closure exceeds tolerance");
    std::ofstream summary(out+"/summary.json");summary<<std::setprecision(15)<<"{\n  \"station\": "<<station<<",\n  \"spacecraft_source\": \""<<source<<"\",\n  \"eop_interpolation\": \"natural cubic spline\",\n  \"eop_units\": \"xp,yp,dX,dY: mas; TAI-UT1: ms\",\n  \"max_closure_ns\": "<<closure*1e9<<",\n  \"max_reference_position_difference_m\": "<<splineError<<",\n  \"rotation_orthogonality\": "<<orthogonality<<"\n}\n";
    std::cout<<"Light-time complete: closure "<<closure*1e9<<" ns, reference difference "<<splineError<<" m\n";
}
#endif
int main(int argc,char**argv){try{
    std::string command=argc>1?argv[1]:"orbit",data="data",out="";std::string trajectorySource="integrated";double hours=20,interval=60,step=60,fitHours=0,countScale=.001;int degree=160,station=34,stride=1,points=180,iterations=12;
    for(int i=2;i<argc;++i){std::string key=argv[i];if(i+1>=argc)throw std::runtime_error("Missing option value");std::string value=argv[++i];if(key=="--data")data=value;else if(key=="--out")out=value;else if(key=="--hours")hours=std::stod(value);else if(key=="--interval")interval=std::stod(value);else if(key=="--degree")degree=std::stoi(value);else if(key=="--step")step=std::stod(value);else if(key=="--trajectory")trajectorySource=value;else if(key=="--count-scale")countScale=std::stod(value);else if(key=="--station")station=std::stoi(value);else if(key=="--stride")stride=std::stoi(value);else if(key=="--points")points=std::stoi(value);else if(key=="--iterations")iterations=std::stoi(value);else if(key=="--fit-hours")fitHours=std::stod(value);else throw std::runtime_error("Unknown option: "+key);}
    if(out.empty())out="results/corrected/"+std::string(command=="fit"?"stage_1_4":command=="doppler"?(trajectorySource=="reference"?"stage_1_3_reference":"stage_1_3_integrated"):command=="lighttime"?"stage_1_2":"stage_1");
    if(command=="selftest"){selftest();return 0;}
    require(command=="orbit"||command=="lighttime"||command=="doppler"||command=="fit","Available commands: orbit, lighttime, doppler, fit, selftest");require(std::isfinite(hours)&&std::isfinite(interval)&&std::isfinite(step)&&std::isfinite(fitHours)&&std::isfinite(countScale)&&countScale>0&&hours>0&&interval>0&&step>0&&degree>=0&&degree<=160&&stride>0&&points>=12&&iterations>=1&&fitHours>=0,"Invalid integration options");
    require(trajectorySource=="integrated"||trajectorySource=="reference","Unknown trajectory source");
    require(command=="selftest"||command=="orbit"||(MESSENGER_STAGE>=2&&command=="lighttime")||(MESSENGER_STAGE>=3&&command=="doppler")||(MESSENGER_STAGE>=4&&command=="fit"),"Command unavailable in this cumulative version");
    require(command!="fit"||trajectorySource=="integrated","Fitting requires the integrated trajectory");
    loadKernels(data);
#if MESSENGER_STAGE >= 2
    if(command=="lighttime"){runLightTime(data,out,jdTDBtoET(2455646.83),hours,interval,station,trajectorySource,step);kclear_c();return 0;}

#endif
#if MESSENGER_STAGE >= 3
    if(command=="doppler"){runDoppler(data,out,stride,countScale,trajectorySource,step);kclear_c();return 0;}

#endif
#if MESSENGER_STAGE >= 4
    if(command=="fit"){runFit(data,out,points,iterations,fitHours,countScale,step);kclear_c();return 0;}

#endif
    require(interval>=step&&std::abs(interval/step-std::round(interval/step))<1e-10,"Output interval must be a multiple of the RK8 step");
    Gravity g(degree);g.load(data+"/HgM009.sha");Dynamics dynamics(g);
    const double epoch=jdTDBtoET(2455646.83);const auto initial=ephemeris(-236,epoch);
    std::filesystem::create_directories(out);std::ofstream csv(out+"/orbit.csv");if(!csv)throw std::runtime_error("Cannot write orbit.csv");
    csv<<"elapsed_s,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s,reference_x_km,reference_y_km,reference_z_km,position_error_m,velocity_error_m_s\n"<<std::setprecision(16);
    Integrator<Dynamics> integ(dynamics,epoch,initial,step);double sum2=0,maxerr=0,finalerr=0;size_t count=0;auto start=std::chrono::steady_clock::now();
    for(double t=0;;t=std::min(t+interval,hours*3600)){const auto state=integ.to(epoch+t),reference=ephemeris(-236,epoch+t);finalerr=1000*norm(pos(state)-pos(reference));maxerr=std::max(maxerr,finalerr);sum2+=finalerr*finalerr;++count;csv<<t;for(double x:state)csv<<','<<x;for(double x:pos(reference))csv<<','<<x;csv<<','<<finalerr<<','<<1000*norm(vel(state)-vel(reference))<<'\n';if(t>=hours*3600)break;}
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream summary(out+"/summary.json");summary<<std::setprecision(15)<<"{\n  \"epoch_jd_tdb\": 2455646.83,\n  \"duration_hours\": "<<hours<<",\n  \"gravity_degree\": "<<degree<<",\n  \"solar_radiation_pressure\": false,\n  \"relativistic_dynamics\": \"EIH beta=gamma=1\",\n  \"integrator\": \"Dormand-Prince RK8 (13 stages)\",\n  \"step_seconds\": "<<step<<",\n  \"GM_mercury_km3_s2\": "<<g.gm<<",\n  \"mercury_orientation\": \"Genova 2019 pole, Stark 2015 rotation formula\",\n  \"position_max_m\": "<<maxerr<<",\n  \"position_rms_m\": "<<std::sqrt(sum2/count)<<",\n  \"position_final_m\": "<<finalerr<<",\n  \"accepted_steps\": "<<integ.stats.accepted<<",\n  \"rejected_steps\": "<<integ.stats.rejected<<",\n  \"runtime_s\": "<<seconds<<"\n}\n";
    std::cout<<std::setprecision(10)<<"Degree "<<degree<<", duration "<<hours<<" h, max error "<<maxerr<<" m, RMS "<<std::sqrt(sum2/count)<<" m, final "<<finalerr<<" m\n"<<"Steps: "<<integ.stats.accepted<<" accepted, "<<integ.stats.rejected<<" rejected; runtime "<<seconds<<" s\n";
    kclear_c();return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
