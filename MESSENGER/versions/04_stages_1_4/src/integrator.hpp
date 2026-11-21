#pragma once
#include "math.hpp"
#include "rk8_table.hpp"
struct IntegrationStats {size_t accepted=0,rejected=0;};
template<class F> class Integrator {
    const F& f_;double time_,step_;State state_;
public:
    IntegrationStats stats;
    Integrator(const F&f,double t,State s,double step=60):f_(f),time_(t),step_(step),state_(s){if(!std::isfinite(step)||step<=0)throw std::runtime_error("Invalid RK8 step");}
    const State& state()const{return state_;}
    State to(double target){
        const double direction=target>=time_?1.:-1.;
        while(direction*(target-time_)>1e-7){
            const double h=direction*std::min(step_,std::abs(target-time_));State k[RK_STAGES]{};
            for(int i=0;i<RK_STAGES;++i){State y=state_;for(int q=0;q<6;++q){double sum=0;for(int j=0;j<i;++j)sum+=RK_A[i][j]*k[j][q];y[q]+=h*sum;}k[i]=f_(time_+RK_C[i]*h,y);}
            for(int q=0;q<6;++q){double sum=0;for(int j=0;j<RK_STAGES;++j)sum+=RK_B[j]*k[j][q];state_[q]+=h*sum;if(!std::isfinite(state_[q]))throw std::runtime_error("Nonfinite RK8 state");}
            time_+=h;++stats.accepted;
        }
        time_=target;return state_;
    }
};
