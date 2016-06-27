//
// Created by dc on 4/8/16.
//

#ifndef RAYTR_PATHTRACER_H
#define RAYTR_PATHTRACER_H

#include "common.h"
#include "scene.h"
#include "image.h"

struct PTParams
{
	int max_depth;
	double p_roulette;
	int emitter_sampling; // 1{Do emitter sampling}
	int n_threads;
	int spp;              // samples per pixel
	int n_rounds;
	int display;
	void LoadFrom (const json &j);
};

class PathTracer
{
	static const int CRD_MAX = 4000, TH_MAX = 64;

	Spectrum **sample_sum;
	int **sample_cnt;

	const Scene *scene_;
	Image *image_;
	int W, H;
	int max_depth_, emitter_sampling_;
	double p_roulette_;
	rand_data rand_datum[TH_MAX];

	int n_threads_, n_samples_, n_rounds_;
	int _probe_x, _probe_y; // For debugging

	Spectrum Trace_ (rand_data *rd, Ray ray, int depth, Spectrum weight, bool last_reflection);
	Spectrum SampleRay (rand_data *rd, int x, int y, int sample_size);

public:
	PathTracer (const Scene &scene, Image &image, const PTParams &p) :
			scene_(&scene), image_(&image),
			max_depth_(p.max_depth), p_roulette_(p.p_roulette),
			emitter_sampling_(p.emitter_sampling),
			n_threads_(p.n_threads), n_samples_(p.spp),
			n_rounds_(p.n_rounds)
	{
		for (int t = 0; t < TH_MAX; ++t)
			rand_init(&rand_datum[t]);

		W = scene.Width();
		H = scene.Height();
		sample_sum = new_matrix<Spectrum>(W, H, Spectrum(0.));
		sample_cnt = new_matrix<int>(W, H, 0);
		if (!emitter_sampling_)
		{
			clog << "Warning: environment light not supported without emitter sampling.\n";
		}
		_probe_x = int(0.37 * W);
		_probe_y = int(0.62 * H);
	}

	~PathTracer ()
	{
		del_matrix<Spectrum>(W, H, sample_sum);
		del_matrix<int>(W, H, sample_cnt);
	}

	Spectrum Trace (rand_data *rd, Ray ray)
	{
		return Trace_(rd, ray, 0, Spectrum(1), false);
	}

	void PT ();
};

#endif //RAYTR_PATHTRACER_H
