//
// Created by dc on 5/17/16.
//

#ifndef RAYTR_BXDF_H
#define RAYTR_BXDF_H

#include "common.h"
#include "random.h"

struct Lambertian: BxDF
{
	virtual double PDF (const Vec3 &I, const Vec3 &N, const Vec3 &O) const override
	{
		// When shading normal is different from geometry normal, it may happen
		// that (I,Ns) > 0. The solution is to treat the light as transmissive
		// instead of reflective. Since the only transmissive BxDF implemented
		// is a Dirac function, we only need to neglect this BxDF.
		if (I.dot(N) > eps)
			return 0.;

		return inv_pi;
	}
	virtual double Sample (rand_data *rd, const Vec3 &, const Vec3 &N, bool, Vec3 &res) const override
	{
		double pdf_;
		CosRandHemisphere(rd, N, res, pdf_);
		return pdf_;
	}

	virtual Type GetType () const override
	{ return BRDF; }
};

struct Blinn_Phong: BxDF
{
	// http://www.thetenthplanet.de/archives/255
private:
	double n, inv_Z;
public:
	Blinn_Phong (double n_): n(n_)
	{
		inv_Z = (n + 2) * (n + 4) / (8 * pi * (n + pow(2.0, -n * 0.5)));
	}
	virtual double PDF (const Vec3 &I, const Vec3 &N, const Vec3 &O) const override
	{
		if (I.dot(N) > eps)
			return 0.;

		Vec3 H = (O - I).normalized();
		return pow(H.dot(N), n) * inv_Z;
	}
	virtual double Sample (rand_data *rd, const Vec3 &I, const Vec3 &N, bool, Vec3 &res) const override
	{
		double costh, pdf;
		Vec3 H, O;
		int ctr = 0;
		do
		{
			// Sample H \sim cos^n<H, N>
			costh = pow(urand(rd), 1.0 / (n + 1));
			H = RandomRotate(rd, sqrt(max(0., 1.0 - sqr(costh))), N);
			O = H * -I.dot(H) * 2 + I;
			++ctr;
		} while (O.dot(N) < eps);
		res = O;
		pdf = pow(costh, n) * (n + 1) / (2 * pi * 4 * O.dot(H));
		return pdf;
	}

	virtual Type GetType () const override
	{ return BRDF; }
};

struct Specular: public BxDF
{
private:
	double refraction;

	static double DielectricReflectance (double costh_i, double n_t_div_i)
	{
		// transmit_div_incident or inside_div_outside
		double sinth_i = sqrt(std::max(0.0, 1. - sqr(costh_i))), sinth_t = sinth_i / n_t_div_i;
		if (sinth_t + eps > 1.) // Total reflectance
			return 1.;

		double costh_t = sqrt(std::max(0.0, 1. - sqr(sinth_t)));
		double ti = n_t_div_i * costh_i, it = costh_t,
				tt = n_t_div_i * costh_t, ii = costh_i;
		double r_par = (ti - it) / (ti + it), r_per = (ii - tt) / (ii + tt);
		return (sqr(r_par) + sqr(r_per)) * 0.5;
	}

public:
	Specular (double refraction_): refraction(refraction_) {}

	virtual double PDF (const Vec3 &I, const Vec3 &N, const Vec3 &O) const override
	{
		throw std::domain_error("Dirac distribution");
	}

	// return 1 if sample is valid, i.e. not rejected due to shading normal effects
	double Sample (rand_data *rd, const Vec3 &I, const Vec3 &N, bool toInside, Vec3 &res) const override
	{
		double costh_i = -I.dot(N);
		double reflectance, refr_ratio;
		Vec3 emergent;

		if (costh_i < -eps)
			return 0.;

		if (refraction > eps)
		{
			refr_ratio = toInside ? refraction : 1.0 / refraction;
			reflectance = DielectricReflectance(costh_i, refr_ratio);
		}
		else
			reflectance = 1.0;

		if (urand(rd) < reflectance + eps)
		{
			// 2.2a specular reflection
			emergent = (I + 2 * costh_i * N).normalized();
		}
		else
		{
			// 2.2b refraction
			emergent = (I + costh_i * N) / refr_ratio;
			emergent = emergent - N * sqrt(std::max(0.0, 1 - emergent.sqnorm()));
		}
		res = emergent;
		return 1.;
	}

	virtual Type GetType () const override
	{ return BxDF_D; }
};

class SpectralSpecular: public BxDF
{
	Specular *bxdfs[SampledSpectrum::N];
public:

	SpectralSpecular (double nD, double disp)
	{
		for (int d = 0; d < SampledSpectrum::N; ++d)
			bxdfs[d] = new Specular((SampledSpectrum::Lambda(d) - 589.3) / (687.0 - 430.8) * disp + nD);
	}

	~SpectralSpecular ()
	{
		for (int d = 0; d < SampledSpectrum::N; ++d)
			delete bxdfs[d];
	}

	virtual double PDF (const Vec3 &I, const Vec3 &N, const Vec3 &O) const override
	{
		throw std::domain_error("SpectralSpecular::PDF called");
	}

	virtual double Sample (rand_data *rd, const Vec3 &I, const Vec3 &N, bool toInside, Vec3 &res) const override
	{
		throw std::domain_error("SpectralSpecular::Sample called");
	}

	// This is a terrible hack
	double Sample (rand_data *rd, const Vec3 &I, const Vec3 &N, bool toInside, SampledSpectrum &Kf, Vec3 &res)
	// @param Kf:	Spectrum accumulated so far (flux / filter). Will be modified by Roulette.
	{
		// Roulette
		double sumPmf[SampledSpectrum::N + 1];
		int item[SampledSpectrum::N + 1], m = 1, j;
		sumPmf[0] = 0;
		for (int d = 0; d < SampledSpectrum::N; ++d)
			if (Kf[d] > eps)
			{
				sumPmf[m] = Kf[d] + (m ? sumPmf[m - 1] : 0);
				item[m] = d;
				m++;
			}

		if (m == 1)
			return 0;

		double psel = urand(rd) * sumPmf[m - 1];
		for (j = m - 1; j > 1 && psel < sumPmf[j - 1]; --j) ;

		// Filter all but j; scale by importance and rejection sampling
		Kf[item[j]] *= (sumPmf[m - 1] / (sumPmf[j] - sumPmf[j - 1]));
		for (int k = 1; k < m; ++k)
			if (k != j) Kf[item[k]] = 0;

		return bxdfs[item[j]]->Sample(rd, I, N, toInside, res);
	}

	virtual Type GetType () const override
	{
		return BxDF_D;
	}
};

#endif //RAYTR_BXDF_H
