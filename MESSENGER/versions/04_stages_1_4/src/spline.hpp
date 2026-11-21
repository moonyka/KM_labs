#pragma once
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
class CubicSpline {
    std::vector<double> x_,y_,second_;
public:
    CubicSpline()=default;
    CubicSpline(std::vector<double>x,std::vector<double>y):x_(std::move(x)),y_(std::move(y)),second_(x_.size()){
        if(x_.size()<3||x_.size()!=y_.size())throw std::runtime_error("Invalid spline data");
        const size_t n=x_.size();std::vector<double> u(n);
        for(size_t i=1;i<n;++i)if(!(x_[i]>x_[i-1]))throw std::runtime_error("Spline epochs must increase");
        for(size_t i=1;i+1<n;++i){const double sigma=(x_[i]-x_[i-1])/(x_[i+1]-x_[i-1]),p=sigma*second_[i-1]+2;second_[i]=(sigma-1)/p;u[i]=(6*((y_[i+1]-y_[i])/(x_[i+1]-x_[i])-(y_[i]-y_[i-1])/(x_[i]-x_[i-1]))/(x_[i+1]-x_[i-1])-sigma*u[i-1])/p;}
        for(size_t k=n-1;k-->0;)second_[k]=second_[k]*second_[k+1]+u[k];
    }
    double operator()(double t)const{
        if(x_.empty()||t<x_.front()||t>x_.back())throw std::runtime_error("Spline epoch outside coverage");
        auto it=std::upper_bound(x_.begin(),x_.end(),t);size_t hi=it==x_.end()?x_.size()-1:size_t(it-x_.begin()),lo=hi-1;
        const double h=x_[hi]-x_[lo],a=(x_[hi]-t)/h,b=(t-x_[lo])/h;
        return a*y_[lo]+b*y_[hi]+((a*a*a-a)*second_[lo]+(b*b*b-b)*second_[hi])*h*h/6;
    }
};
