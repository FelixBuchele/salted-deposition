#include "sdbd/orientation.hpp"
#include <cmath>
#include <stdexcept>

namespace sdbd {
Quaternion Quaternion::normalized() const {
    const auto length=std::hypot(std::hypot(w,x),std::hypot(y,z));
    if (!std::isfinite(length) || length==0) throw std::invalid_argument("Invalid orientation quaternion");
    return {w/length,x/length,y/length,z/length};
}
Vec3 Quaternion::apply(Vec3 value) const {
    const Vec3 axis(x,y,z);
    const auto cross=2*v(axis,value);
    return value+w*cross+v(axis,cross);
}
Quaternion operator*(const Quaternion& a,const Quaternion& b) {
    return Quaternion{a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z,
        a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y,
        a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x,
        a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w}.normalized();
}
Quaternion Quaternion::from_axis(Vec3 axis,long double cosine,long double sine) {
    const auto half=std::atan2(sine,cosine)/2;
    const auto length=norm(axis);
    if (!std::isfinite(length) || length==0) throw std::invalid_argument("Invalid rotation axis");
    const auto factor=std::sin(half)/length;
    return Quaternion{std::cos(half),factor*axis.x(),factor*axis.y(),factor*axis.z()}.normalized();
}
Quaternion Quaternion::from_basis(const std::array<Vec3,3>& c) {
    long double m[3][3]={{c[0].x(),c[1].x(),c[2].x()},
                         {c[0].y(),c[1].y(),c[2].y()},
                         {c[0].z(),c[1].z(),c[2].z()}};
    std::array<long double,4> q{};
    const auto trace=m[0][0]+m[1][1]+m[2][2];
    if (trace>0) {
        const auto t=2*std::sqrt(1+trace);
        q={t/4,(m[2][1]-m[1][2])/t,(m[0][2]-m[2][0])/t,(m[1][0]-m[0][1])/t};
    } else {
        int i=m[1][1]>m[0][0] ? 1 : 0;
        if (m[2][2]>m[i][i]) i=2;
        const int j=(i+1)%3,k=(i+2)%3;
        const auto t=2*std::sqrt(1+m[i][i]-m[j][j]-m[k][k]);
        q[0]=(m[k][j]-m[j][k])/t;q[i+1]=t/4;
        q[j+1]=(m[j][i]+m[i][j])/t;q[k+1]=(m[k][i]+m[i][k])/t;
    }
    return Quaternion{q[0],q[1],q[2],q[3]}.normalized();
}
}
