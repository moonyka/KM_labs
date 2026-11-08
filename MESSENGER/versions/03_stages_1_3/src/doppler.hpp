#pragma once
#include "lighttime.hpp"
#include <cstdint>
#include <iostream>
#include <iomanip>

inline constexpr double OBS_JD_BASE=2455643.0;
inline double taiOffset(double jd1,double jd2){return ((jd1-OBS_JD_BASE)+jd2)*86400.;}
struct Observation {double tai,count,reference,observed;int tx,rx;std::array<double,6> extra;bool eligible=true;};
inline std::vector<Observation> readObservations(const std::string&path,double countScale=0.001){
    if(!std::isfinite(countScale)||countScale<=0)throw std::runtime_error("Count scale must be positive");
    std::ifstream f(path);size_t n;f>>n;if(!f||n==0)throw std::runtime_error("Invalid observation header");
    std::vector<Observation> result;std::string line;std::getline(f,line);
    for(size_t i=0;i<n;++i){if(!std::getline(f,line))throw std::runtime_error("Truncated observation file");std::istringstream r(line);double j1,j2,compression,high,low,oh,ol;int type;Observation o{};
        r>>j1>>j2>>type>>compression>>high>>low>>o.rx>>o.tx>>oh>>ol;if(!r)throw std::runtime_error("Invalid observation row "+std::to_string(i+2));
        // The teacher's text export uses milliseconds; ignore columns after obs_frac.
        if(type!=12)continue;o.tai=taiOffset(j1,j2);o.count=compression*countScale;o.reference=(high*16777216.+low)*.001;o.observed=oh+ol*1e-9;
        if(!std::isfinite(o.tai)||!std::isfinite(o.count)||!std::isfinite(o.reference)||!std::isfinite(o.observed)||o.count<=0)throw std::runtime_error("Invalid type 12 observation");result.push_back(o);
    }
    std::sort(result.begin(),result.end(),[](const auto&a,const auto&b){return a.tai<b.tai;});if(result.empty())throw std::runtime_error("No type 12 observations");
    std::vector<double> values;for(const auto&o:result)values.push_back(o.observed);std::sort(values.begin(),values.end());const double median=values[values.size()/2];for(auto&v:values)v=std::abs(v-median);std::sort(values.begin(),values.end());const double cutoff=100*std::max(values[values.size()/2],.001);
    for(auto&o:result)o.eligible=std::abs(o.observed-median)<=cutoff;
    return result;
}
class Ramps {
    struct Segment {double begin,end,rate,frequency;};std::map<int,std::vector<Segment>> segments_;
public:
    explicit Ramps(const std::string&path){std::ifstream f(path);if(!f)throw std::runtime_error("Missing ramp table");std::string line;
        while(std::getline(f,line)){if(line.empty()||line[0]=='d')continue;std::istringstream r(line);double a,b,ri,rf,fg,fi,ff,c,d;int station,duplicate;r>>a>>b>>station>>ri>>rf>>fg>>fi>>ff>>duplicate>>c>>d;if(!r)throw std::runtime_error("Invalid ramp row");Segment s{taiOffset(a,b),taiOffset(c,d),ri+rf*1e-9,fg*1e9+std::fmod(fi,1e9)+ff*1e-9};if(station!=duplicate||s.end<=s.begin||fg==0)throw std::runtime_error("Invalid ramp interval or non-sky-level frequency");segments_[station].push_back(s);}
        for(auto&[id,v]:segments_)std::sort(v.begin(),v.end(),[](const auto&a,const auto&b){return a.begin<b.begin;});
    }
    double phase(int station,double start,double duration)const{
        const auto it=segments_.find(station);if(it==segments_.end())throw std::runtime_error("Missing ramp station");double total=0,remaining=duration,t=start;
        while(remaining>1e-10){const auto&v=it->second;auto s=std::upper_bound(v.begin(),v.end(),t+1e-8,[](double x,const Segment&b){return x<b.begin;});if(s==v.begin())throw std::runtime_error("Transmission precedes ramp coverage");--s;if(t<s->begin-1e-7||t>=s->end+1e-7)throw std::runtime_error("Gap in ramp coverage");const double dt=std::min(remaining,s->end-t);if(dt<=0)throw std::runtime_error("Invalid ramp boundary");const double f=s->frequency+s->rate*(t-s->begin);total+=f*dt+.5*s->rate*dt*dt;t+=dt;remaining-=dt;}
        return total;
    }
};
inline double etMinusTai(double et,int station,const Stations&stations){
    const double ttMinusTdb=ephemeris(1000000001,et,1000000000)[0];
    const Vec earthVelocity=vel(ephemeris(399,et,0)),site=stations.geocentric(station,et);
    return 32.184-ttMinusTdb+dot(earthVelocity,site)/(LIGHT_SPEED*LIGHT_SPEED);
}
inline double taiOffsetToET(double offset){const double tt=(OBS_JD_BASE-2451545.)*86400.+offset+32.184;double et=tt;for(int i=0;i<3;++i)et=tt-ephemeris(1000000001,et,1000000000)[0];return et;}
inline double stationTAItoET(double offset,int station,const Stations&stations){const double tai=(OBS_JD_BASE-2451545.)*86400.+offset;double et=taiOffsetToET(offset);for(int i=0;i<3;++i)et=tai+etMinusTai(et,station,stations);return et;}
struct DopplerResult {double computed,residual,roundtrip;};
inline DopplerResult computeDoppler(const Observation&o,const Ramps&ramps,const LightTime&solver,const Stations&stations){
    const double a=o.tai-o.count/2,b=o.tai+o.count/2,eta=stationTAItoET(a,o.rx,stations),etb=stationTAItoET(b,o.rx,stations);
    const auto pa=solver.solve(eta,o.tx,o.rx),pb=solver.solve(etb,o.tx,o.rx);const double taua=pa.up_s+pa.down_s,taub=pb.up_s+pb.down_s;
    const double da=etMinusTai(eta,o.rx,stations),db=etMinusTai(etb,o.rx,stations),dea=etMinusTai(pa.transmitET(),o.tx,stations),deb=etMinusTai(pb.transmitET(),o.tx,stations);
    // Keep durations and small time differences separate from large J2000 epochs.
    const double txStart=a-taua+(da-dea),txDuration=o.count-(taub-taua)+(db-da)-(deb-dea);
    const double received=ramps.phase(o.tx,txStart,txDuration)/o.count,ratio=880./749.;const double computed=ratio*(o.reference-received);
    return {computed,o.observed-computed,.5*(taua+taub)};
}
inline void runDoppler(const std::string&data,const std::string&out,int stride,double countScale=0.001){
    const auto observations=readObservations(data+"/obs_4.txt",countScale);const Ramps ramps(data+"/ramps_4.txt");const Stations stations(data);
    const double begin=taiOffsetToET(observations.front().tai-observations.front().count),end=taiOffsetToET(observations.back().tai+observations.back().count);
    SpacecraftSpline trajectory(begin,end);LightTime solver(stations,[&](double et){return trajectory(et);});
    std::filesystem::create_directories(out);std::ofstream f(out+"/doppler.csv");f<<std::setprecision(16)<<"tai_offset_s,station,count_s,reference_hz,observed_hz,computed_hz,residual_hz,roundtrip_s,eligible_for_fit\n";double sum=0,sum2=0,maxabs=0,rawSum2=0;size_t used=0,skipped=0,flagged=0;
    std::ofstream rejected(out+"/rejected.csv");rejected<<"tai_offset_s,station,reason\n";
    for(size_t i=0;i<observations.size();i+=stride){const auto&o=observations[i];try{const auto r=computeDoppler(o,ramps,solver,stations);f<<o.tai<<','<<o.rx<<','<<o.count<<','<<o.reference<<','<<o.observed<<','<<r.computed<<','<<r.residual<<','<<r.roundtrip<<','<<o.eligible<<'\n';rawSum2+=r.residual*r.residual;if(o.eligible){sum+=r.residual;sum2+=r.residual*r.residual;maxabs=std::max(maxabs,std::abs(r.residual));++used;}else{++flagged;rejected<<std::setprecision(16)<<o.tai<<','<<o.rx<<",extreme_observable_100_MAD\n";}}catch(const std::runtime_error&e){++skipped;rejected<<std::setprecision(16)<<o.tai<<','<<o.rx<<",\""<<e.what()<<"\"\n";}}
    if(!used)throw std::runtime_error("No computable Doppler observations; inspect rejected.csv");
    std::ofstream summary(out+"/summary.json");summary<<std::setprecision(15)<<"{\n  \"spacecraft_source\": \"msgr_final.bsp (reference-orbit validation)\",\n  \"count_epoch\": \"TAI midpoint (ODF convention)\",\n  \"count_scale_seconds\": "<<countScale<<",\n  \"eop_epoch_scale\": \"TAI as specified by teacher\",\n  \"tt_tdb_source\": \"de442t.bsp time ephemeris\",\n  \"reference_convention\": \"Moyer f_REF = (880/749) * decoded REF_FREQ\",\n  \"turnaround_ratio\": \"880/749\",\n  \"total_type12\": "<<observations.size()<<",\n  \"stride\": "<<stride<<",\n  \"used\": "<<used<<",\n  \"flagged_extreme_observations\": "<<flagged<<",\n  \"outlier_rule\": \"abs(observed - median) > 100 * MAD; original rows retained in doppler.csv\",\n  \"rejected_computations\": "<<skipped<<",\n  \"raw_residual_rms_hz\": "<<std::sqrt(rawSum2/(used+flagged))<<",\n  \"residual_mean_hz\": "<<sum/used<<",\n  \"residual_rms_hz\": "<<std::sqrt(sum2/used)<<",\n  \"residual_max_abs_hz\": "<<maxabs<<"\n}\n";
    std::cout<<"Doppler: used "<<used<<", rejected "<<skipped<<", mean "<<sum/used<<" Hz, RMS "<<std::sqrt(sum2/used)<<" Hz, max "<<maxabs<<" Hz\n";
}
