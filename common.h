//
// common.h: external headers & common tools
// (Vector, Matrix, Ray, Maybe, new/del-matrix)
// Created by dc on 4/7/16.
//

#ifndef RAYTR_MATH_H
#define RAYTR_MATH_H

// #define DIRECT_ONLY
// #define SPECTRAL

#include "stdafx.h"

// {{{ Math

const double eps = 1e-7, inf = 1e30;
const double pi = acos(-1), inv_pi = 1.0 / pi;

template <typename T>
inline T sqr (T x) { return x * x; }

struct Vec3
{
	double x, y, z;

	Vec3 (double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_)
	{ }

	double operator[] (int id) const
	{ return (&x)[id]; }

	double dot (const Vec3 &b) const
	{ return x * b.x + y * b.y + z * b.z; }

	double sqnorm () const
	{ return sqr(x) + sqr(y) + sqr(z); }

	double norm () const
	{ return sqrt(sqnorm()); }

	Vec3 operator+ (const Vec3 &b) const
	{ return Vec3(x + b.x, y + b.y, z + b.z); }

	Vec3 operator- (const Vec3 &b) const
	{ return Vec3(x - b.x, y - b.y, z - b.z); }

	Vec3 operator- () const
	{ return Vec3(-x, -y, -z); }

	Vec3 operator* (double l) const
	{ return Vec3(l * x, l * y, l * z); }

	Vec3 operator/ (double l) const
	{ return Vec3(x / l, y / l, z / l); }

	Vec3 &operator+= (const Vec3 &b)
	{
		x += b.x;
		y += b.y;
		z += b.z;
		return *this;
	}

	Vec3 normalized ()
	{
		double l = sqnorm();
		assert(l > sqr(eps));
		return *this / sqrt(l);
	}

	Vec3 cross (const Vec3 &b) const
	{
		return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
	}

	static Vec3 Parse (const json &j)
	{
		assert(j.is_array() && j.size() == 3);
		return Vec3(j[0].get<double>(), j[1].get<double>(), j[2].get<double>());
	}
};

using CPVec3 = const Vec3 *;

inline std::ostream& operator<< (std::ostream& out, const Vec3 &b)
{
	return out << '(' << b.x << ", " << b.y << ", " << b.z << ')';
}

inline Vec3 operator* (double l, const Vec3 &b)
{
	return Vec3(l * b.x, l * b.y, l * b.z);
}

inline void GetPerpVecs (const Vec3 &N, Vec3 &A, Vec3 &B)
{
	if (fabs(N.y) + fabs(N.x) < 1e-2)
	{
		// numeric stability
		A = Vec3(0, N.z, -N.y);
		B = Vec3(-(sqr(N.y) + sqr(N.z)), N.x * N.y, N.x * N.z);
	}
	else
	{
		A = Vec3(N.y, -N.x, 0);
		B = Vec3(N.z * N.x, N.z * N.y, -(sqr(N.x) + sqr(N.y)));
	}
	A = A.normalized();
	B = B.normalized();
}

struct Ray
{
	Vec3 S, D;

	Ray () {}

	Ray (Vec3 S_, Vec3 D_) : S(S_), D(D_)
	{ }
};

struct VecTf
{
	Eigen::Matrix4d m;

	VecTf (const Eigen::Matrix4d &m_) : m(m_)
	{ }

	VecTf (const json &json_)
	{
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				m(i, j) = json_.at(i).at(j).get<double>();
	}

	static VecTf Eye ()
	{ return VecTf(Eigen::Matrix4d::Identity()); }

	Vec3 operator() (const Vec3 &v) const
	{
		Eigen::Vector4d vv;
		vv << v.x, v.y, v.z, 1.;
		vv = m * vv;
		return Vec3(vv(0), vv(1), vv(2)) / vv(3);
	}

	Vec3 TVec (const Vec3 &v) const
	{
		Eigen::Vector4d vv;
		vv << v.x, v.y, v.z, 0.;
		vv = m * vv;
		return Vec3(vv(0), vv(1), vv(2));
	}

	VecTf Inverse () const
	{ return VecTf(m.inverse()); }
};

// }}} Math

template <typename T>
struct Maybe
{
	bool valid;
	T v;
	static Maybe<T> Some (const T &b) { return Maybe<T>{true, b}; }
	static Maybe<T> None () { return Maybe<T>{false, T()}; }
	operator bool () const { return valid; }
};

template <typename T> T** new_matrix (int n, int m, const T &v)
{
	T **ret = new T*[n];
	for (int i = 0; i < n; ++i)
	{
		ret[i] = new T[m];
		std::fill(ret[i], ret[i] + m, v);
	}
	return ret;
}
template <typename T> void del_matrix (int n, int m, T **v)
{
	for (int i = 0; i < n; ++i)
		delete[] v[i];

	delete[] v;
}

template <typename T>
struct Averager
{
	T sum;
	int cnt;

	Averager (T init) : sum(init), cnt(0)
	{ }

	void Add (T v)
	{
		sum += v;
		cnt++;
	}

	T Get ()
	{ return sum / (double) cnt; }

	void Reset (const T &init)
	{
		sum = init;
		cnt = 0;
	}
};

// Handles resource that the owner doesn't need to know the details of.
// What is the preferred practice for this?
class Destroyable
{
public:
	virtual ~Destroyable () {}
};

inline double Interp (double xs, double xe, double ys, double ye, double xp)
{
	return ys + (ye - ys) * (xp - xs) / (xe - xs);
}

static double Average (const double x[], const double y[], int n, double x0, double x1)
{
	// Corner cases
	if (x1 < x[0]) return y[0];
	if (x0 > x[n - 1]) return y[n - 1];

	double s = 0;
	// CC: [?, 0) and [n-1, ?)
	if (x0 < x[0]) s += y[0] * (x[0] - x0);
	if (x1 > x[n - 1]) s += y[n - 1] * (x1 - x[n - 1]);

	for (int d = 0; d + 1 < n; ++d)
	{
		// Add contribution of [d,d+1)
		double xs = x[d], xe = x[d + 1], ys = y[d], ye = y[d + 1];
		double vs = max(x0, xs), ve = min(x1, xe);
		if (vs > ve) continue;
		s += (Interp(xs, xe, ys, ye, vs) + Interp(xs, xe, ys, ye, ve)) * (ve - vs) / 2;
	}
	return s / (x1 - x0);
}

#endif //RAYTR_MATH_H
