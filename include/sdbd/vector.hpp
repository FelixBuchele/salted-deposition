// Adapted from Nikola Topic's SDBD source (2015); vector code credits Thomas Schwager.
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include <cmath>
#include <iostream>

namespace sdbd {
// Three Cartesian long-double components. Coordinates use the caller's length
// unit; vector products retain the original numerical expression order.
class Vec3 {
    // Read x,y,z in stream order; propagate the stream's extraction state.
    friend std::istream& operator>>(std::istream& input, Vec3& value) {
        return input >> value.x_ >> value.y_ >> value.z_;
    }
    // Write space-separated x,y,z using the caller's precision and locale.
    friend std::ostream& operator<<(std::ostream& output, const Vec3& value) {
        return output << value.x_ << ' ' << value.y_ << ' ' << value.z_;
    }
    // Return componentwise addition without changing either operand.
    friend Vec3 operator+(const Vec3& a, const Vec3& b) {
        Vec3 result(a); result+=b; return result;
    }
    // Return componentwise subtraction without changing either operand.
    friend Vec3 operator-(const Vec3& a, const Vec3& b) {
        Vec3 result(a); result-=b; return result;
    }
    // Return value multiplied componentwise by scalar c.
    friend Vec3 operator*(long double c, const Vec3& value) {
        Vec3 result(value); result*=c; return result;
    }
    // Return scalar multiplication with the same ordering as c*value.
    friend Vec3 operator*(const Vec3& value, long double c) { return c*value; }
    // Return squared Euclidean norm, in squared component units.
    friend long double norm2(const Vec3& value) {
        return value.x_*value.x_+value.y_*value.y_+value.z_*value.z_;
    }
    // Return Euclidean norm in component units; caller ensures numerical range.
    friend long double norm(const Vec3& value) { return std::sqrt(norm2(value)); }
    // Return Cartesian scalar product, retaining left-to-right addition order.
    friend long double s(const Vec3& a, const Vec3& b) {
        return a.x_*b.x_+a.y_*b.y_+a.z_*b.z_;
    }
    // Return right-handed Cartesian cross product a x b.
    friend Vec3 v(const Vec3& a, const Vec3& b) {
        Vec3 result;
        result.x_=a.y_*b.z_-a.z_*b.y_;
        result.y_=a.z_*b.x_-a.x_*b.z_;
        result.z_=a.x_*b.y_-a.y_*b.x_;
        return result;
    }
    // Return Rodrigues rotation of r about unit axis n through rp; cc and ss
    // are the angle cosine/sine. The caller supplies a normalized valid axis.
    friend Vec3 rot(const Vec3& r,const Vec3& rp,long double cc,long double ss,const Vec3& n) {
        Vec3 rc=r-rp;
        return (1-cc)*s(n,rc)*n+cc*rc+ss*v(n,rc)+rp;
    }
public:
    // Construct Cartesian components; an omitted component is zero.
    constexpr explicit Vec3(long double x=0,long double y=0,long double z=0):x_(x),y_(y),z_(z) {}
    // Return mutable x component; units follow the represented quantity.
    long double& x() { return x_; }
    // Return x component by value.
    long double x() const { return x_; }
    // Return mutable y component; units follow the represented quantity.
    long double& y() { return y_; }
    // Return y component by value.
    long double y() const { return y_; }
    // Return mutable z component; units follow the represented quantity.
    long double& z() { return z_; }
    // Return z component by value.
    long double z() const { return z_; }
    // Add value componentwise in place and return this vector.
    Vec3& operator+=(const Vec3& value) {
        x_+=value.x_; y_+=value.y_; z_+=value.z_; return *this;
    }
    // Subtract value componentwise in place and return this vector.
    Vec3& operator-=(const Vec3& value) {
        x_-=value.x_; y_-=value.y_; z_-=value.z_; return *this;
    }
    // Multiply components by c in place and return this vector.
    Vec3& operator*=(long double c) { x_*=c; y_*=c; z_*=c; return *this; }
private:
    long double x_,y_,z_;
};

inline constexpr Vec3 null(0,0,0);
inline constexpr Vec3 ex(1,0,0);
inline constexpr Vec3 ey(0,1,0);
inline constexpr Vec3 ez(0,0,1);
}
