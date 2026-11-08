#pragma once
#include "math.hpp"
#include <vector>
// Dormand-Prince 5(4), with separate position and velocity absolute tolerances.
struct IntegrationStats {size_t accepted=0,rejected=0;};
template<class F> class Integrator {
    const F& f_;double time_,step_;State state_;double tolerance_;
public:
    IntegrationStats stats;
    Integrator(const F&f,double t,State s,double tolerance=1e-8):f_(f),time_(t),step_(10),state_(s),tolerance_(tolerance){}
    const State& state()const{return state_;}
    State to(double target){
        if(target<time_)throw std::runtime_error("Integrator requires nondecreasing epochs");
        static const double c[7]={0,1./5,3./10,4./5,8./9,1,1};
        static const double a[7][7]={{0},{1./5},{3./40,9./40},{44./45,-56./15,32./9},{19372./6561,-25360./2187,64448./6561,-212./729},{9017./3168,-355./33,46732./5247,49./176,-5103./18656},{35./384,0,500./1113,125./192,-2187./6784,11./84}};
        static const double bh[7]={35./384,0,500./1113,125./192,-2187./6784,11./84,0};
        static const double bl[7]={5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};
        while(time_<target){double h=std::min(step_,target-time_);if(h<1e-7){time_=target;break;}State k[7];
            for(int i=0;i<7;++i){State y=state_;for(int j=0;j<i;++j)for(int q=0;q<6;++q)y[q]+=h*a[i][j]*k[j][q];k[i]=f_(time_+c[i]*h,y);}
            State next=state_;double error=0;
            for(int q=0;q<6;++q){double e=0;for(int j=0;j<7;++j){next[q]+=h*bh[j]*k[j][q];e+=h*(bh[j]-bl[j])*k[j][q];}const double scale=tolerance_*(q<3?1.:.001)+1e-13*std::max(std::abs(state_[q]),std::abs(next[q]));error=std::max(error,std::abs(e)/scale);}
            if(!std::isfinite(error))throw std::runtime_error("Nonfinite integration error");
            if(error<=1){state_=next;time_+=h;++stats.accepted;}else ++stats.rejected;
            const double factor=error==0?5:std::clamp(.9*std::pow(error,-.2),.1,5.);
            step_=std::min(120.,h*factor);if(step_<1e-5)throw std::runtime_error("Integration step underflow");
        }return state_;
    }
};
