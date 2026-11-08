#pragma once
#include "math.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Real, fully (4-pi) normalized solid harmonics, without Condon-Shortley phase.
// Cartesian recurrence avoids the longitude singularity on the polar axis.
class Gravity {
    int degree_;
    std::vector<double> c_,s_,a_,b_;
    mutable std::vector<Jet> v_,w_;
    static size_t ix(int n,int m){return size_t(n)*(n+1)/2+m;}
public:
    double gm=22031.863566,radius=2440;
    explicit Gravity(int n=160):degree_(n),c_(ix(n+1,0)),s_(c_.size()),a_(c_.size()),b_(c_.size()),v_(c_.size()),w_(c_.size()){
        c_[0]=1;
        for(int m=0;m<=n;++m)for(int l=m+2;l<=n;++l){
            const double ll=double(l)*l,mm=double(m)*m;
            a_[ix(l,m)]=std::sqrt((4*ll-1)/(ll-mm));
            b_[ix(l,m)]=std::sqrt((2.*l+1)*((l-1.)*(l-1)-mm)/((2.*l-3)*(ll-mm)));
        }
    }
    int degree()const{return degree_;}
    void set(int n,int m,double c,double s=0){c_.at(ix(n,m))=c;s_.at(ix(n,m))=s;}
    void load(const std::string& path){
        std::ifstream f(path);if(!f)throw std::runtime_error("Cannot open gravity file: "+path);
        std::string line;std::getline(f,line);std::replace(line.begin(),line.end(),',',' ');std::istringstream h(line);
        double unused;int n,m,normalization;h>>radius>>gm>>unused>>n>>m>>normalization;
        if(!h||n<degree_||m<degree_||normalization!=1)throw std::runtime_error("Unsupported gravity header/degree/normalization");
        size_t loaded=0;
        while(std::getline(f,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream r(line);int l,k;double c,s;if(!(r>>l>>k>>c>>s))continue;if(l<0||k<0||k>l)throw std::runtime_error("Invalid gravity coefficient");if(l<=degree_){set(l,k,c,s);++loaded;}}
        // SHA files omit degree zero: the central term is explicitly retained.
        c_[0]=1;
        if(loaded+1<c_.size())throw std::runtime_error("Incomplete gravity coefficients");
    }
    Jet potential(const Vec& p)const{
        const Jet x=Jet::variable(p[0],0),y=Jet::variable(p[1],1),z=Jet::variable(p[2],2);
        const Jet rr=x*x+y*y+z*z;
        if(rr.v<1)throw std::runtime_error("Invalid gravity position");
        const Jet r=sqrt(rr), X=Jet(radius)*x/rr,Y=Jet(radius)*y/rr,Z=Jet(radius)*z/rr,Q=Jet(radius*radius)/rr;
        v_[0]=Jet(radius)/r;w_[0]=Jet(0);
        Jet sum=v_[0];
        for(int m=0;m<=degree_;++m){
            if(m>0){const auto prev=ix(m-1,m-1),k=ix(m,m);const double d=m==1?std::sqrt(3.):std::sqrt((2.*m+1)/(2.*m));
                v_[k]=Jet(d)*(X*v_[prev]-Y*w_[prev]);w_[k]=Jet(d)*(X*w_[prev]+Y*v_[prev]);}
            for(int n=m;n<=degree_;++n){const auto k=ix(n,m);
                if(n==m+1){v_[k]=Jet(std::sqrt(2.*m+3))*Z*v_[ix(m,m)];w_[k]=Jet(std::sqrt(2.*m+3))*Z*w_[ix(m,m)];}
                if(n>=m+2){v_[k]=Jet(a_[k])*Z*v_[ix(n-1,m)]-Jet(b_[k])*Q*v_[ix(n-2,m)];w_[k]=Jet(a_[k])*Z*w_[ix(n-1,m)]-Jet(b_[k])*Q*w_[ix(n-2,m)];}
                if(n>0)sum=sum+Jet(c_[k])*v_[k]+Jet(s_[k])*w_[k];
            }
        }
        return Jet(gm/radius)*sum;
    }
    Vec acceleration(const Vec&p)const{return potential(p).d;}
};
