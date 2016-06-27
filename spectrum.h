//
// Created by dc on 5/20/16.
//

#ifndef RAYTR_SPECTRUM_H
#define RAYTR_SPECTRUM_H

#include "common.h"


struct RGBSpectrum
{
	static const int N = 3;
	double r, g, b; // [0, 1]
	RGBSpectrum (double r_, double g_, double b_) : r(r_), g(g_), b(b_)
	{ }

	RGBSpectrum (double s = 0) : r(s), g(s), b(s)
	{ }

	double operator[] (int id) const
	{ return (&r)[id]; }
	double& operator[] (int id)
	{ return (&r)[id]; }

	RGBSpectrum operator+ (const RGBSpectrum &u) const
	{ return RGBSpectrum(r + u.r, g + u.g, b + u.b); }

	RGBSpectrum operator- (const RGBSpectrum &u) const
	{ return RGBSpectrum(r - u.r, g - u.g, b - u.b); }

	RGBSpectrum operator* (const RGBSpectrum &u) const
	{ return RGBSpectrum(r * u.r, g * u.g, b * u.b); }

	RGBSpectrum operator/ (double u) const
	{ return RGBSpectrum(r / u, g / u, b / u); }

	RGBSpectrum &operator*= (const RGBSpectrum &u)
	{
		r *= u.r;
		g *= u.g;
		b *= u.b;
		return *this;
	}

	RGBSpectrum &operator+= (const RGBSpectrum &u)
	{
		r += u.r;
		g += u.g;
		b += u.b;
		return *this;
	}

	bool nonzero () const
	{ return r > eps || g > eps || b > eps; }

	RGBSpectrum ToRGB () const { return *this; }

	double y () const { return 0.21267 * r + 0.71516 * g + 0.07217 * b; }

	static RGBSpectrum Parse (const json &j)
	{
		if (j.count("v"))
			return RGBSpectrum(J_AT_F(j, "v"), J_AT_F(j, "v"), J_AT_F(j, "v"));
		else
		{
			assert(j.is_array() && j.size() == 3);
			return RGBSpectrum(j[0].get<double>(), j[1].get<double>(), j[2].get<double>());
		}
	}
};

inline RGBSpectrum operator* (double l, const RGBSpectrum b)
{
	return RGBSpectrum(l * b.r, l * b.g, l * b.b);
}

inline std::ostream& operator<< (std::ostream& out, const RGBSpectrum &b)
{
	return out << '(' << b.r << ", " << b.g << ", " << b.b << ')';
}

struct SampledSpectrum
{
	static const int N = 28; // 400 - 700
	static SampledSpectrum X, Y, Z;
	static double yint;
	typedef Eigen::Matrix<float, N, 1> Vecf;
	Vecf vs;

	static constexpr double DLambda () { return (700 - 400) / double(N); }
	static constexpr double Lambda (int n)
	{
		return 400. + 300 * ((double)n / N) + DLambda() / 2;
	}
	static constexpr double Lambda0 (int n)
	{
		return 400. + 300 * ((double)n / N);
	}

	static void Init ();

	SampledSpectrum () { }

	SampledSpectrum (double i): vs(Vecf::Constant(i)) {}
	SampledSpectrum (const Vecf &v): vs(v) {}
	SampledSpectrum (const SampledSpectrum &r): vs(r.vs) {}

	float& operator[] (int v) { return vs[v]; }
	float operator[] (int v) const { return vs[v]; }

	SampledSpectrum operator+ (const SampledSpectrum &u) const
	{
		return SampledSpectrum(vs + u.vs);
	}

	SampledSpectrum operator- (const SampledSpectrum &u) const
	{
		return SampledSpectrum(vs - u.vs);
	}

	SampledSpectrum operator* (const SampledSpectrum &u) const
	{
		return SampledSpectrum(vs.cwiseProduct(u.vs));
	}

	SampledSpectrum operator* (double u) const
	{
		return SampledSpectrum(vs * u);
	}

	SampledSpectrum operator/ (double u) const
	{
		return (*this) * (1. / u);
	}

	SampledSpectrum& operator*= (const SampledSpectrum &u)
	{
		vs = vs.cwiseProduct(u.vs);
		return *this;
	}

	SampledSpectrum& operator+= (const SampledSpectrum &u)
	{
		vs += u.vs;
		return *this;
	}

	bool nonzero () const
	{
		for (int d = 0; d < N; ++d)
			if (vs[d] > eps)
				return true;

		return false;
	}

	RGBSpectrum ToRGB () const
	{
		double xyz[3] = {0, 0, 0};
		RGBSpectrum rgb;
		for (int d = 0; d < N; ++d)
		{
			xyz[0] += vs[d] * X[d] * DLambda();
			xyz[1] += vs[d] * Y[d] * DLambda();
			xyz[2] += vs[d] * Z[d] * DLambda();
		}
		xyz[0] /= yint;
		xyz[1] /= yint;
		xyz[2] /= yint;
		rgb[0] =  3.240479*xyz[0] - 1.537150*xyz[1] - 0.498535*xyz[2];
		rgb[1] = -0.969256*xyz[0] + 1.875991*xyz[1] + 0.041556*xyz[2];
		rgb[2] =  0.055648*xyz[0] - 0.204043*xyz[1] + 1.057311*xyz[2];

		return RGBSpectrum(rgb[0], rgb[1], rgb[2]);
	}

	double y () const
	{
		return vs.cwiseProduct(Y.vs).sum() * DLambda() / yint;
	}

	static SampledSpectrum FromRGB (RGBSpectrum rgb, bool reflective = true);

	static SampledSpectrum Parse (const json &j, bool reflective = true)
	{
		if (j.count("v"))
		{
			return SampledSpectrum(j.at("v").get<double>());
		}
		assert(j.is_array());
		if (j.size() == 3)
		{
			// RGB
			return FromRGB(RGBSpectrum::Parse(j), reflective);
		}
		SampledSpectrum sp;
		const int input_N = 30;
		double rs[input_N], x[input_N];
		for (int d = 0; d < input_N; ++d)
		{
			rs[d] = j[d].get<double>();
			x[d] = 400 + 300.0 * d / input_N + 5;
		}
		for (int d = 0; d < N; ++d)
			sp[d] = Average(x, rs, input_N, Lambda(d), Lambda(d + 1));

		return sp;
	}
};

inline SampledSpectrum operator* (double l, const SampledSpectrum &b)
{
	return b * l;
}

inline std::ostream& operator<< (std::ostream& out, const SampledSpectrum &b)
{
	out << "[";
	for (int d = 0; d < SampledSpectrum::N; ++d)
		out << b[d] << ", ";
	return out << "]";
}


#ifndef SPECTRAL
typedef RGBSpectrum Spectrum;
#else
typedef SampledSpectrum Spectrum;
#endif

#endif //RAYTR_SPECTRUM_H
