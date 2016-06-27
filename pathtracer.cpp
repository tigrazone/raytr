//
// Created by dc on 4/8/16.
//

#include "pathtracer.h"
#include "shade.h"
#include "bxdfs.h"

const double eps_mv = 5e-3; // eps_mv ** 2 > eps

void PathTracer::PT ()
{
//	SampleRay(rand_datum, 245, 600 - 443, 10);

	std::atomic<int> counter(0);
	int cur_round;
	auto th_disp = std::thread([&]() {
		while (cur_round < n_rounds_)
		{
			cerr << "\r" << "Round " << cur_round << ": " << counter << "/" <<
			(uint64_t)W * H;
			timespec tim, tim2;
			tim.tv_nsec = int(1e8);
			tim.tv_sec = 0;
			nanosleep(&tim, &tim2);
		}
	});
	for (cur_round = 0; cur_round < n_rounds_; ++cur_round)
	{
		vector<std::thread> threads(n_threads_);
		vector<uint32_t> samples((uint64_t)W * H);

		for (uint32_t x = 0, pos = 0; x < W; ++x)
		{
			for (uint32_t y = 0; y < H; ++y)
				samples[pos++] = (x << 16u) | y;
		}

		std::random_shuffle(samples.begin(), samples.end());
		using iter_t = decltype(samples.begin());
		counter = 0;

		auto sample_thread = [this, &counter](rand_data *rd, iter_t s, iter_t e)
		{
			int ctr = 0;
			for (iter_t i = s; i < e; ++i)
			{
				int x = i[0] >> 16u, y = i[0] & 0xffff;
				Spectrum color = SampleRay(rd, x, y, n_samples_);
				image_->Set(x, y, color);
				if (++ctr % 200 == 0)
				{
					counter += ctr;
					ctr = 0;
				}
			}
			counter += ctr;
		};

		for (int t = 0; t < n_threads_; ++t)
		{
			size_t div = samples.size() / n_threads_;
			size_t s = div * t, e = div * (t + 1);
			if (t + 1 == n_threads_)
				e = samples.size();

			iter_t is = samples.begin() + s, ie = samples.begin() + e;

			threads[t] = std::thread([&sample_thread, t, is, ie, this](){
				sample_thread(&rand_datum[t], is, ie);
			});
		}
		for (auto &th: threads)
			th.join();

		clog << "C(" << _probe_x << ", " << _probe_y << ") = "
		     << std::ios::fixed << std::setprecision(4)
		     << sample_sum[_probe_x][_probe_y] / (cur_round + 1) << endl;
		image_->Dump(cur_round);
	}
	th_disp.join();
}

Spectrum PathTracer::SampleRay (rand_data *rd, int x, int y, int sample_size)
{
	const Camera* camera = scene_->GetCamera();
	Spectrum &sum = sample_sum[x][y];
    Spectrum cur(0);
	for (int t = 0; t < sample_size; ++t)
	{
		double jtrx = urand(rd), jtry = urand(rd);
		Ray ray;
		camera->Project(rd, x + jtrx, y + jtry, ray);
		cur = cur + Trace(rd, ray);
	}
    sum += cur.operator*(1.0 / sample_size);
	sample_cnt[x][y] ++;
	return sum.operator*(1.0 / sample_cnt[x][y]);
};

Spectrum PathTracer::Trace_ (rand_data *rd, Ray ray, int depth, Spectrum weight, bool last_specular)
{
	if (depth > this->max_depth_)
	{
		return Spectrum(0);
	}

	auto tmp = scene_->GetCollision(ray);
	if (!tmp)
		return Spectrum(0);

	Spectrum ret(0);
	auto collision = tmp.v;
	auto newpos = collision.pos;
	auto surface = collision.obj->GetSurfaceProperty();

	// 1. LOCAL EMISSION
	// - include iff first hit or (last ray is specular (no emitter sample was drawn))
	if (!emitter_sampling_ || (depth == 0 || last_specular))
		if (surface->emission.nonzero())
			ret = ret + weight * surface->emission;

	// 2. CHOOSE: diffusion - specular - absorption
	auto bxdf_i = surface->SampleBxDF(rd);
	auto bxdf = bxdf_i->bxdf;
	Spectrum K = bxdf_i->GetColor(collision.obj, collision.pos);

	if (bxdf == nullptr) // emitter
		goto PT_END;

	if (! bxdf->IsDirac())
	{
		// 2.1.1: direct light. here we sample all lights.
		if (emitter_sampling_)
		{
			Spectrum colorLS(0.), colorBS(0.);
			// a. env light
			colorLS += K * scene_->GetEnvLight();

			// b. others
			const auto &lights = scene_->GetLights();
			for (auto l: lights)
			{
				auto light = l.first;
				double light_pmf = l.second;
				if (light == collision.obj)
					continue;

				Vec3 lsPoint, lsDir;
				lsPoint = light->SampleSurface(rd, newpos);
				lsDir = (lsPoint - newpos).normalized();

				// a. Check if sampled light can reach collision.pos
				if (lsDir.dot(collision.norm) < 0)
					// collision.obj separates observer and light
					continue;

				tmp = scene_->GetCollision(Ray(newpos, lsDir));
				if (tmp.v.obj != light)
					// sheltered by other objects
					continue;

				assert((tmp.v.pos - lsPoint).sqnorm() < 1e-2 ||
					   (tmp.v.pos - newpos).cross(lsPoint - newpos).sqnorm() < 1e-4);

				/* b. Shade
				 * I = \int L_e(point) * bxdf_pdf * cos(theta) d\omega
				 * \hat{I} = L_e(point) * bxdf_pdf / area_light_pdf * cos(theta)
				 * ret += light_pmf * K * I
				 */
				double light_invpdf = light->SAInvPDF(newpos);
				double bxdf_pdf = bxdf->PDF(ray.D, collision.norm, lsDir);

				Spectrum I = light->GetSurfaceProperty()->emission;
				I *= light_pmf * bxdf_pdf * light_invpdf * lsDir.dot(collision.norm);
				I *= K;
				colorLS = colorLS + I;
			}

			// Check MIS someday.
			ret = ret + weight * colorLS;
		}
#ifndef DIRECT_ONLY
		// 2.1.2: recurse
		double dir_pdf;
		Vec3 nD;
		dir_pdf = bxdf->Sample(rd, ray.D, collision.norm, 0, nD);

		Spectrum nW = weight
			/ dir_pdf * nD.dot(collision.norm)
			* bxdf->PDF(ray.D, collision.norm, nD)
			* K;
		double p_roulette = depth > 5 ? 
			                min(p_roulette_, nW.y()) :
			                p_roulette_;
		if (urand(rd) < p_roulette)
		{
			Ray nRay(newpos + eps_mv * nD, nD);
			nW *= 1.0 / p_roulette;
			ret = ret + Trace_(rd, nRay, depth + 1, nW, false);
		}
#endif
	}
	else // Dirac BxDF
	{
		// 2.2.0 Get reflectance
		Spectrum specular(0);
		Vec3 nD;
		double valid;
#ifndef SPECTRAL
		valid = bxdf->Sample(rd, ray.D, collision.norm, collision.inside, nD);
#else
		auto bxdfs = dynamic_cast<SpectralSpecular*>(bxdf);
		valid = bxdfs->Sample(rd, ray.D, collision.norm, collision.inside, weight, nD);
#endif
		Spectrum kW = weight * K;
		double p_roulette = depth > 5 ? 
			                min(p_roulette_, kW.y()) :
			                p_roulette_;
		if (valid > eps && urand(rd) < p_roulette)
		{
			Ray nRay(newpos + eps_mv * nD, nD);
			specular = Trace_(rd, nRay, depth + 1, (1.0 / p_roulette) * kW, true);
		}
		ret = ret + specular;
	}

PT_END:
	assert(!(-1 * ret).nonzero());
	// 2.3: absorption
	return ret;
}

void PTParams::LoadFrom (const json &j)
{
#define L(d, X, v) X=J_AT_##d##_NZ(j, #X, v)
	L(I, max_depth, 11);
	L(F, p_roulette, 0.98);
	L(I, emitter_sampling, 1);
	L(I, n_threads, 4);
	L(I, spp, 10);
	L(I, n_rounds, 10);
#ifndef NO_CV
	L(I, display, 1);
#else
	L(I, display, 0);
#endif
}
