#include "lighttime.hpp"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv)
{
 loadKernels(argv[1]); Stations stations(argv[1]);
 LightTime light(stations,[](Epoch t){return pos(ephemeris(-236,t,0));});
 std::ifstream file(std::string(argv[1])+"/lt_4.txt");std::string header;std::getline(file,header);
 long double jd,tx; double down,up;int station,row=0;
 std::cout<<std::setprecision(16)<<"row,down_difference_ns,up_difference_ns\n";
 while(file>>jd>>down>>up>>tx>>station){
  if(row%1000==0){
   const long double day=std::floor(jd);
   Epoch t(static_cast<double>((day-2451545.L)*86400.L),static_cast<double>((jd-day)*86400.L));
   auto model=light.solve(t,station,station);
   std::cout<<row<<','<<(model.down_s-down)*1e9<<','<<(model.up_s-up)*1e9<<'\n';
  }
  ++row;
 }
}
