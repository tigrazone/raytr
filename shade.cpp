//
// Created by dc on 5/16/16.
//

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"

#include "shade.h"
#include "bxdfs.h"

Texture::Texture (const string &path)
{
#ifndef SPECTRAL
	int comp_got;
	uint8_t *text = stbi_load(path.c_str(), &w, &h, &comp_got, 4);
	if (comp_got < 3 || comp_got > 4 || !text)
		throw std::logic_error("Texture not loaded: " + path);

	data = new_matrix<Spectrum>(w, h, Spectrum(0.));
	for (int x = 0; x < w; ++x)
	{
		for (int y = 0; y < h; ++y)
		{
			int p = x * h + y;
			data[x][y] = Spectrum(text[p * comp_got + 0] / 255.0,
							   text[p * comp_got + 1] / 255.0,
							   text[p * comp_got + 2] / 255.0);
		}
	}
#else
	throw std::domain_error("texture not supported in spectral mode");
#endif
}

inline double Intensity (Spectrum &s)
{
	double i(0);
	for (int d = 0; d < Spectrum::N; ++d)
		i = max(i, (double)s[d]);

	return i;
}

ShadingProperty::ShadingProperty (Spectrum diffusion_, Spectrum emission_,
								  Spectrum specular_,
								  double phong_exp,
								  double refraction_,
								  const string &mKd_path,
								  const string &mKs_path,
								  const string &bump_path,
								  const std::map<string, string> *kwargs)
{
	emission = emission_;
	n_bxdfs = 0;
	if (diffusion_.nonzero() || mKd_path != "")
	{
		bxdfs[n_bxdfs++] = BxDFItem{
				new Lambertian(),
				diffusion_, Texture::Create(mKd_path),
				Intensity(diffusion_)};
	}

	if (specular_.nonzero() || mKs_path != "")
	{
		auto mKs_map = Texture::Create(mKs_path);
		if (phong_exp < eps || phong_exp > 1e4)
		{
#ifdef SPECTRAL
			double nD = kwargs && kwargs->count("nD") ? 
				atof(kwargs->at("nD").c_str()) : 1.;
			double disp = kwargs && kwargs->count("dispD") ?
				atof(kwargs->at("dispD").c_str()) : 0.;
			bxdfs[n_bxdfs++] = BxDFItem{
					new SpectralSpecular(nD, disp),
					specular_, nullptr,
					Intensity(specular_)
			};
#else
			bxdfs[n_bxdfs++] = BxDFItem{
					new Specular(refraction_),
					specular_, mKs_map,
					Intensity(specular_)
			};
#endif
		}
		else
		{
			bxdfs[n_bxdfs++] = BxDFItem{
					new Blinn_Phong(phong_exp),
					specular_, mKs_map,
					Intensity(specular_)
			};
		}
	}

	if (emission_.nonzero())
	{
		bxdfs[n_bxdfs++] = BxDFItem{
				nullptr, emission_, Texture::Create(""),
				Intensity(emission_)
		};
	}

	assert(n_bxdfs > 0);

	double sum = 0;
	for (int t = 0; t < n_bxdfs; ++t)
		sum += bxdfs[t].pmf;

	for (int t = 0; t < n_bxdfs; ++t)
		bxdfs[t].pmf /= sum + eps;

// map_bump = Texture::Create(bump_path); TODO
}

PSurfaceProperty ShadingProperty::Create (const json &j)
{
	// In json only p_diffusion and p_specular is defined.
	Spectrum color = Spectrum::Parse(j.at("diffusion"));
	Spectrum diffusion = color * J_AT_F(j, "p_diffusion");
	Spectrum specular = color * J_AT_F(j, "specular");
	Spectrum emission = j.count("emission") ? Spectrum::Parse(j.at("emission")) : Spectrum(0);
	std::map<string, string> kw;
	for (json::const_iterator it = j.begin(); it != j.end(); ++it)
	{
		if (it.value().is_string())
			kw[it.key()] = it.value();
	}
	auto ret = make_shared<ShadingProperty>(
			diffusion,
			emission,
			specular,
			J_AT_F_NZ(j, "phong_exp", 0.0),
			J_AT_F_NZ(j, "refraction", 0.0),
			j.count("mKd") ? j["mKd"].get<string>() : "",
			j.count("mKs") ? j["mKs"].get<string>() : "",
			j.count("bump_texname") ? j["bump_texname"].get<string>() : "",
			&kw
	);
	return ret;
}


void LightProperty::SampleDir (rand_data *rd, Vec3 &pos, Vec3 &N, Vec3 &dir, double &pos_invpdf, double &dir_pdf) const
{
	prim_->SampleSurface(rd, pos, N);
	pos_invpdf = prim_->Area();
	// TODO: PT doesn't support these hacks yet.
	if (light_exp_ == 1 && !unidirected_light_)
		CosRandHemisphere(rd, N, dir, dir_pdf);
	else if (light_exp_ != 1)
	{
		dir = RandomRotate(rd, sqrt(1.0 - pow(urand(rd), 2.0 / (1 + light_exp_))), N);
		dir_pdf = 1.;
	}
	else
	{
		N = light_i_;
		double costh = pow(urand(rd), 0.001),
				sinth = sqrt(max(0., 1.0 - sqr(costh))) * (urand(rd) > 0.5 ? 1 : -1);
		dir = light_i_ * costh + light_j_ * sinth;
		dir_pdf = 1.;
	}
}

