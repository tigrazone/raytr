//
// shade.h: texture and BxDFs
// Created by dc on 4/13/16.
//

#ifndef RAYTR_SURFACE_H
#define RAYTR_SURFACE_H

#include "common.h"
#include "random.h"
#include "objekt.h"
#include "spectrum.h"

struct Texture
{
	Spectrum **data;
	int w, h;

	Texture (const string &path);

	Spectrum at (Vec3 P) const
	{
		double u = P.x, v = P.y;
		u = u - floor(u);
		v = v - floor(v);
		int px = int(u * (w - 1)), py = int(v * (h - 1));
		double rx = u * (w - 1) - px, ry = v * (h - 1) - py;
		auto ret = rx * ry * data[px + 1][py + 1] +
				   (1 - rx) * ry * data[px][py + 1] +
				   rx * (1 - ry) * data[px + 1][py] +
				   (1 - rx) * (1 - ry) * data[px][py];
		return ret;
	}

	static shared_ptr<Texture> Create (const string &path)
	{
		if (path == "")
			return shared_ptr<Texture>(nullptr);

		return make_shared<Texture>(path);
	}
};

struct BxDF
{
public:
	enum Type { BRDF, BxDF_D };
	// I points to surface; O points out of surface
	virtual double PDF (const Vec3 &I, const Vec3 &N, const Vec3 &O) const = 0;
	virtual double Sample (rand_data *rd, const Vec3 &I, const Vec3 &N, bool toInside, Vec3 &res) const = 0;
	virtual Type GetType () const = 0;
	inline bool IsDirac () const { return GetType() == BxDF_D; }
};

struct BxDFItem
{
	BxDF *bxdf;
	Spectrum K;
	shared_ptr<Texture> Kmap;
	double pmf;

	Spectrum GetColor (const Primitive *obj, const Vec3 &Pworld) const
	{
		if (Kmap == nullptr)
			return K / pmf;

		return Kmap->at(obj->ToUV(Pworld));
	}
};

struct ShadingProperty
{
 	BxDFItem bxdfs[4];
	int n_bxdfs;
	Spectrum emission;

	// NOTE: if emission_ > 0, direct lighting expects Kd, Ks etc. should be small.
	// Otherwise result will be inconsistent.
	ShadingProperty (Spectrum diffusion_, Spectrum emission_,
					 Spectrum specular_,
					 double phong_exp,
					 double refraction_,
					 const string &mKd_path = "",
					 const string &mKs_path = "",
					 const string &bump_path = "",
					 const std::map<string, string> *kwargs = nullptr);

	const BxDFItem *SampleBxDF (rand_data *rd) const
	{
		double p = urand(rd);
		for (int t = 0; t < n_bxdfs; ++t)
		{
			if (p < bxdfs[t].pmf || t + 1 == n_bxdfs)
				return &bxdfs[t];
			else
				p -= bxdfs[t].pmf;
		}
	}
/*
	double SumBxDF (BxDF::Type type, const Vec3 &I, const Vec3 &N, const Vec3 &O)
	{
		double f = 0, pmf = 0;
		for (int t = 0; t < n_bxdfs; ++t)
			if (bxdfs[t].bxdf && bxdfs[t].bxdf->GetType() == type)
			{
				f += bxdfs[t].bxdf->PDF(I, N, O);
				pmf += bxdfs[t].pmf;
			}

		return f / pmf;
	}*/

	static PSurfaceProperty Create (const json &j);

};

struct LightProperty
{
	const Primitive *prim_;
	Vec3 light_i_, light_j_;
	bool unidirected_light_;
	int light_exp_;

	LightProperty (): prim_(nullptr), unidirected_light_(false), light_exp_(1) {}

	void SampleDir (rand_data *rd, Vec3 &pos, Vec3 &N, Vec3 &dir, double &pos_invpdf, double &dir_pdf) const;
	inline bool HasPDF () const
	{
		return !unidirected_light_ && light_exp_ == 1;
	}
};

#endif //RAYTR_SURFACE_H
