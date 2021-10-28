#pragma once
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

#define     GLH_ZERO                double(0.0)
#define     GLH_EPSILON          double(10e-6)
#define		GLH_EPSILON_2		double(10e-12)
#define     equivalent(a,b)             (((a < b + GLH_EPSILON) &&\
                                                      (a > b - GLH_EPSILON)) ? true : false)
inline double lerp(double a, double b, float t)
{
	return a + t * (b - a);
}

inline double fmax(double a, double b) {
	return (a > b) ? a : b;
}

inline double fmin(double a, double b) {
	return (a < b) ? a : b;
}

inline bool isEqual(double a, double b, double tol = GLH_EPSILON)
{
	return fabs(a - b) < tol;
}

#ifndef M_PI
#define M_PI 3.14159f
#endif

class vec3f {
public:
	union {
		struct {
			double x, y, z;
		};
		struct {
			double v[3];
		};
	};

	inline vec3f()
	{
		x = 0; y = 0; z = 0;
	}

	inline vec3f(const vec3f &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	inline vec3f(const double *v)
	{
		x = v[0];
		y = v[1];
		z = v[2];
	}

	inline vec3f(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	inline double operator [] (int i) const { return v[i]; }
	inline double &operator [] (int i) { return v[i]; }

	inline vec3f &operator += (const vec3f &v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	inline vec3f &operator = (const vec3f &v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	inline vec3f &operator -= (const vec3f &v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	inline vec3f &operator *= (double t) {
		x *= t;
		y *= t;
		z *= t;
		return *this;
	}

	inline vec3f &operator /= (double t) {
		x /= t;
		y /= t;
		z /= t;
		return *this;
	}

	inline void negate() {
		x = -x;
		y = -y;
		z = -z;
	}

	inline vec3f operator - () const {
		return vec3f(-x, -y, -z);
	}

	inline vec3f operator+ (const vec3f &v) const
	{
		return vec3f(x + v.x, y + v.y, z + v.z);
	}

	inline vec3f operator- (const vec3f &v) const
	{
		return vec3f(x - v.x, y - v.y, z - v.z);
	}

	inline vec3f operator *(double t) const
	{
		return vec3f(x*t, y*t, z*t);
	}

	inline vec3f operator /(double t) const
	{
		return vec3f(x / t, y / t, z / t);
	}

	// cross product
	inline const vec3f cross(const vec3f &vec) const
	{
		return vec3f(y*vec.z - z * vec.y, z*vec.x - x * vec.z, x*vec.y - y * vec.x);
	}

	inline double dot(const vec3f &vec) const {
		return x * vec.x + y * vec.y + z * vec.z;
	}

	inline void normalize()
	{
		double sum = x * x + y * y + z * z;
		if (sum > GLH_EPSILON_2) {
			double base = double(1.0 / sqrt(sum));
			x *= base;
			y *= base;
			z *= base;
		}
	}

	inline double length() const {
		return double(sqrt(x*x + y * y + z * z));
	}

	inline vec3f getUnit() const {
		return (*this) / length();
	}

	inline bool isUnit() const {
		return isEqual(squareLength(), 1.f);
	}

	//! max(|x|,|y|,|z|)
	inline double infinityNorm() const
	{
		return fmax(fmax(fabs(x), fabs(y)), fabs(z));
	}

	inline vec3f & set_value(const double &vx, const double &vy, const double &vz)
	{
		x = vx; y = vy; z = vz; return *this;
	}

	inline bool equal_abs(const vec3f &other) {
		return x == other.x && y == other.y && z == other.z;
	}

	inline double squareLength() const {
		return x * x + y * y + z * z;
	}

	static vec3f zero() {
		return vec3f(0.f, 0.f, 0.f);
	}

	//! Named constructor: retrieve vector for nth axis
	static vec3f axis(int n) {
		assert(n < 3);
		switch (n) {
		case 0: {
			return xAxis();
		}
		case 1: {
			return yAxis();
		}
		case 2: {
			return zAxis();
		}
		}
		return vec3f();
	}

	//! Named constructor: retrieve vector for x axis
	static vec3f xAxis() { return vec3f(1.f, 0.f, 0.f); }
	//! Named constructor: retrieve vector for y axis
	static vec3f yAxis() { return vec3f(0.f, 1.f, 0.f); }
	//! Named constructor: retrieve vector for z axis
	static vec3f zAxis() { return vec3f(0.f, 0.f, 1.f); }

};

inline vec3f operator * (double t, const vec3f &v) {
	return vec3f(v.x*t, v.y*t, v.z*t);
}

inline vec3f interp(const vec3f &a, const vec3f &b, double t)
{
	return a * (1 - t) + b * t;
}

inline vec3f vinterp(const vec3f &a, const vec3f &b, double t)
{
	return a * t + b * (1 - t);
}

inline vec3f interp(const vec3f &a, const vec3f &b, const vec3f &c, double u, double v, double w)
{
	return a * u + b * v + c * w;
}

inline double clamp(double f, double a, double b)
{
	return fmax(a, fmin(f, b));
}

inline double vdistance(const vec3f &a, const vec3f &b)
{
	return (a - b).length();
}


inline std::ostream& operator<<(std::ostream&os, const vec3f &v) {
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
	return os;
}

inline void
vmin(vec3f &a, const vec3f &b)
{
	a.set_value(
		fmin(a[0], b[0]),
		fmin(a[1], b[1]),
		fmin(a[2], b[2]));
}

inline void
vmax(vec3f &a, const vec3f &b)
{
	a.set_value(
		fmax(a[0], b[0]),
		fmax(a[1], b[1]),
		fmax(a[2], b[2]));
}

inline vec3f lerp(const vec3f &a, const vec3f &b, float t)
{
	return a + t * (b - a);
}


inline double fmax(double a, double b, double c)
{
	double t = a;
	if (b > t) t = b;
	if (c > t) t = c;
	return t;
}

inline double fmin(double a, double b, double c)
{
	double t = a;
	if (b < t) t = b;
	if (c < t) t = c;
	return t;
}

inline int project3(const vec3f &ax,
	const vec3f &p1, const vec3f &p2, const vec3f &p3)
{
	double P1 = ax.dot(p1);
	double P2 = ax.dot(p2);
	double P3 = ax.dot(p3);

	double mx1 = fmax(P1, P2, P3);
	double mn1 = fmin(P1, P2, P3);

	if (mn1 > 0) return 0;
	if (0 > mx1) return 0;
	return 1;
}

inline int project6(vec3f &ax,
	vec3f &p1, vec3f &p2, vec3f &p3,
	vec3f &q1, vec3f &q2, vec3f &q3)
{
	double P1 = ax.dot(p1);
	double P2 = ax.dot(p2);
	double P3 = ax.dot(p3);
	double Q1 = ax.dot(q1);
	double Q2 = ax.dot(q2);
	double Q3 = ax.dot(q3);

	double mx1 = fmax(P1, P2, P3);
	double mn1 = fmin(P1, P2, P3);
	double mx2 = fmax(Q1, Q2, Q3);
	double mn2 = fmin(Q1, Q2, Q3);

	if (mn1 > mx2) return 0;
	if (mn2 > mx1) return 0;
	return 1;
}

