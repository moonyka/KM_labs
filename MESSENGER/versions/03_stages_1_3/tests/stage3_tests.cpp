#include "../src/doppler.hpp"
#include <chrono>

void expect(bool value,const char* message){if(!value)throw std::runtime_error(message);}
struct TempFile {
    std::filesystem::path path;
    explicit TempFile(const char* kind):path(std::filesystem::temp_directory_path()/(std::string("messenger_")+kind+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".txt")){}
    ~TempFile(){std::error_code ec;std::filesystem::remove(path,ec);}
};
int main(){try{
    TempFile obs("observations_");
    {std::ofstream f(obs.path);f<<"2\n2455643 0 12 500 0 1000 34 34 -12 250000000 ignored columns\n2455643 0 11 500 0 1000 34 34 0 0\n";}
    const auto observations=readObservations(obs.path.string());
    expect(observations.size()==1,"Type 12 filter failed");
    expect(observations[0].count==.5,"Teacher milliseconds must give 0.5 seconds for 500");
    expect(observations[0].reference==1&&observations[0].observed==-11.75,"Packed observables failed");
    expect(readObservations(obs.path.string(),.01)[0].count==5,"Explicit ODF scale failed");
    bool invalid=false;try{readObservations(obs.path.string(),0);}catch(const std::runtime_error&){invalid=true;}expect(invalid,"Zero count scale accepted");
    TempFile ramp("ramps_");
    {std::ofstream f(ramp.path);f<<std::setprecision(17)<<"d 34 2\n2455643 0 34 2 0 7 1000000005 0 34 2455643 "<<10./86400<<"\n2455643 "<<10./86400<<" 34 -1 0 7 25 0 34 2455643 "<<20./86400<<"\n";}
    const Ramps ramps(ramp.path.string());
    expect(std::abs(ramps.phase(34,0,.5)-3500000002.75)<1e-5,"Linear ramp integral or frequency modulus failed");
    expect(std::abs(ramps.phase(34,9,3)-21000000072.)<1e-4,"Cross-boundary phase integral failed");
    expect(std::abs(ramps.phase(34,10,2)-14000000048.)<1e-4,"Exact ramp boundary failed");
    bool gap=false;try{ramps.phase(34,19,2);}catch(const std::runtime_error&){gap=true;}expect(gap,"Missing ramp coverage accepted");
    bool missing=false;try{ramps.phase(99,0,1);}catch(const std::runtime_error&){missing=true;}expect(missing,"Missing station accepted");
    std::cout<<"PASS: teacher milliseconds, observation fields, ignored columns, type filter, ramp frequency modulus, analytic phase and boundary/gap handling\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
