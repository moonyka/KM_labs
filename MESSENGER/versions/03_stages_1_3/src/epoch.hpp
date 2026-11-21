#pragma once
#include <cmath>
struct Epoch {
    double high=0,low=0;
    Epoch()=default;
    Epoch(double value):high(value){}
    Epoch(double a,double b){high=a+b;const double z=high-a;low=(a-(high-z))+(b-z);}
    Epoch shifted(double seconds)const{return Epoch(high,low+seconds);}
};
inline double difference(const Epoch&a,const Epoch&b){return (a.high-b.high)+(a.low-b.low);}
