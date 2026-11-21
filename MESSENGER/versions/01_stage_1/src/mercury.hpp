#pragma once
#include "math.hpp"
inline constexpr double DEG_TO_RAD=3.14159265358979323846/180.;
inline void mercuryRotation(double et,double rotation[3][3]){
    const double days=et/86400.,centuries=days/36525.;
    const double alpha=(281.0082-.032808*centuries)*DEG_TO_RAD,delta=(61.4164-.0048464*centuries)*DEG_TO_RAD;
    const double amplitude[5]={1.,-.105231253725171,-.0103439039634868,-.0023786959435392,-.000535306724141346};
    double w=329.5469+6.1385025*days;
    for(int k=1;k<=5;++k)w+=(40./3600.)*amplitude[k-1]*std::sin(k*4.09233445*(days+42.71182)*DEG_TO_RAD);
    w=std::remainder(w,360.)*DEG_TO_RAD;
    const double ca=std::cos(alpha),sa=std::sin(alpha),cd=std::cos(delta),sd=std::sin(delta),cw=std::cos(w),sw=std::sin(w);
    const Vec equatorX{-sa,ca,0},equatorY{-sd*ca,-sd*sa,cd},pole{cd*ca,cd*sa,sd};
    const Vec x=cw*equatorX+sw*equatorY,y=(-sw)*equatorX+cw*equatorY;
    for(int i=0;i<3;++i){rotation[0][i]=x[i];rotation[1][i]=y[i];rotation[2][i]=pole[i];}
}
