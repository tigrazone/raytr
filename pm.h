//
// Created by dc on 4/12/16.
//

#ifndef RAYTR_PM_H
#define RAYTR_PM_H

#include <mutex>
#include <vector>
#include <utility>
#include "common.h"
#include "objekt.h"
#include "scene.h"
#include "image.h"

struct PMPixel
{
	int scr_x, scr_y;
	Averager<Spectrum> color;
	PMPixel (): color(0)
	{
	}
};

struct PMRay
{
	Vec3 P, N, D;
	BxDF *bxdf;
	Spectrum weight; // [eye, obj]
	double sqrad, cnt_photons;
	Spectrum tau;
	PMPixel *pixel;
	PMRay (PMPixel *pixel_): pixel(pixel_)
	{
		tau = Spectrum(0.);
		sqrad = cnt_photons = 0.;
	}
};

struct PhotonInfo
{
	Vec3 dir;
	Spectrum flux; // (obj, light]
	double dis; // distance traveled
};

template <typename Elem>
struct ParallelVector
{
	std::vector<Elem> vec;
	std::mutex mutex_;
	void Clear ()
	{
		std::lock_guard<std::mutex> lock(mutex_);
		vec.clear();
	}
	void Merge (const std::vector<Elem> &part)
	{
		std::lock_guard<std::mutex> lock(mutex_);
		vec.insert(vec.end(), part.begin(), part.end());
	}
};

struct PMParams
{
	int max_depth, ray_max_depth;
	int start_round;
	double R0; // initial photon radius
	double p_caustic;
	double p_roulette;
	double alpha;
	int n_threads, n_rays, photons_per_round, n_rounds;
	int display;
	void LoadFrom (const json &j);
};

class PM
{
	static const int TH_MAX = 48;
	static const int
			PM_CAUSTIC = 1,  // Caustic path
			PM_NON_CAUSTIC = 2, // Any path except caustic/direct light
			PM_DIRECT_LIGHT = 4; // Direct light path

	const double PM_ALPHA;
	const Scene *scene_;
	Image *image_;
	ParallelVector<std::pair<Vec3, PhotonInfo>> photons_;
	ParallelVector<PMRay> rays_;
	std::vector<std::thread> threads_;
	double total_photons_;
	Spectrum **photon_map_;
	PMPixel **pixels_;
	int **n_ray_samples_, n_ray_sample_base_;
	rand_data rand_datum[TH_MAX];

	int max_depth_, width_, height_, start_round_, ray_max_depth_;
	double p_caustic_, p_roulette_, R0_;
	int n_threads_, n_rays_, photons_per_round_, n_rounds_;

	int _probe_x, _probe_y; // for debugging

	void _UpdatePMImage (const vector<pair<Vec3, PhotonInfo>> &photons);
	void _PreSampleRays (rand_data *rd, int n_threads, int n_rays);

	void EmitRays (int n_threads);
	void EmitPhotons (int n_threads, int n_photons);
	bool CastRay (rand_data *rd, const Ray &ray, PMRay &res);
	void CastPhoton (rand_data *rd, Vec3 src, PhotonInfo photon, vector<pair<Vec3, PhotonInfo>> &pool,
						 int castMode);

public:
	PM (const Scene &scene, Image &image,
		const PMParams &p):
		scene_(&scene), image_(&image),
		max_depth_(p.max_depth), ray_max_depth_(p.ray_max_depth),
		start_round_(p.start_round), R0_(p.R0),
		p_caustic_(p.p_caustic), p_roulette_(p.p_roulette),
		n_threads_(p.n_threads), n_rays_(p.n_rays),
		photons_per_round_(p.photons_per_round),
		n_rounds_(p.n_rounds), PM_ALPHA(p.alpha)
	{
		if (!ray_max_depth_) ray_max_depth_ = max_depth_;

		for (int t = 0; t < TH_MAX; ++t)
			rand_init(&rand_datum[t]);

		total_photons_ = 0.;
		width_ = scene.Width();
		height_ = scene.Height();

		photon_map_ = new Spectrum*[width_];
		for (int w = 0; w < width_; ++w)
		{
			photon_map_[w] = new Spectrum[height_];
			std::fill(photon_map_[w], photon_map_[w] + height_, Spectrum(0.));
		}

		pixels_ = new_matrix<PMPixel>(width_, height_, PMPixel());
		for (int x = 0; x < width_; ++x)
			for (int y = 0; y < height_; ++y)
			{
				pixels_[x][y].scr_x = x;
				pixels_[x][y].scr_y = y;
			}

		_probe_x = int(0.37 * width_);
		_probe_y = int(0.62 * height_);
	}
	void PPM ();
	void PPPM ();
};


#endif //RAYTR_PM_H
