#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include <stdexcept>
using Vec = std::array<double,3>;
using State = std::array<double,6>;
inline Vec operator+(Vec a,const Vec& b){for(int i=0;i<3;++i)a[i]+=b[i];return a;}
inline Vec operator-(Vec a,const Vec& b){for(int i=0;i<3;++i)a[i]-=b[i];return a;}
inline Vec operator*(double s,Vec a){for(auto& x:a)x*=s;return a;}
inline double dot(const Vec&a,const Vec&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline double norm(const Vec&a){return std::sqrt(dot(a,a));}
inline Vec pos(const State&s){return {s[0],s[1],s[2]};}
inline Vec vel(const State&s){return {s[3],s[4],s[5]};}
inline Vec mv(const double a[3][3],const Vec&v){Vec r{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)r[i]+=a[i][j]*v[j];return r;}
inline Vec mtv(const double a[3][3],const Vec&v){Vec r{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)r[i]+=a[j][i]*v[j];return r;}
// First derivatives of a scalar with respect to Cartesian x, y, z.
struct Jet {double v=0;Vec d{};Jet()=default;Jet(double x):v(x){};static Jet variable(double x,int k){Jet a(x);a.d[k]=1;return a;}};
inline Jet operator+(const Jet&a,const Jet&b){Jet r(a.v+b.v);r.d=a.d+b.d;return r;}
inline Jet operator-(const Jet&a,const Jet&b){Jet r(a.v-b.v);r.d=a.d-b.d;return r;}
inline Jet operator*(const Jet&a,const Jet&b){Jet r(a.v*b.v);r.d=b.v*a.d+a.v*b.d;return r;}
inline Jet operator/(const Jet&a,const Jet&b){Jet r(a.v/b.v);r.d=(1/b.v)*(a.d-(a.v/b.v)*b.d);return r;}
inline Jet sqrt(const Jet&a){Jet r(std::sqrt(a.v));r.d=(0.5/r.v)*a.d;return r;}
