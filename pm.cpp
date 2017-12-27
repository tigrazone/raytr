//
// Created by dc on 4/12/16.
//

#include "pm.h"
#include "photon_tree.h"
#include "shade.h"
#include "bxdfs.h"

using namespace std;

static const double eps_mv = 5e-3;

inline bool _GetPixel (const Scene *scene_, Vec3 pnt, int &x, int &y)
{
	Maybe<Vec3> r = scene_->GetCamera()->Rasterize(pnt);
	x = int(r.v.x);
	y = int(r.v.y);
	return r.valid;
}

void PM::_UpdatePMImage (const vector<pair<Vec3, PhotonInfo>> &photons)
{
	static std::mutex mut;
	std::lock_guard<std::mutex> guard(mut);
	for (auto ph: photons)
	{
		int x, y;
		if (_GetPixel(scene_, ph.first, x, y))
		{
			photon_map_[x][y] += ph.second.flux;
			auto c = photon_map_[x][y];
			image_->Set(width_ + x, y, c);
		}
	}
}

void PM::EmitPhotons (int n_threads, int n_photons)
{
	total_photons_ += n_photons;

	// Emit photon for a range of light samples
	auto emitPhotons = [&](rand_data *rd, int cur_n_photons){
		vector<pair<Vec3, PhotonInfo>> pool;
		int cnt = 0, merge_every = 20000;
		for (int i = 0; i < cur_n_photons; ++i)
		{
			// Choose light with roulette
			const Primitive *c_light;
			double light_pmf;
			c_light = scene_->SampleLight(rd, light_pmf);

			// Calc scaled flux
			Vec3 P_lit, N_lit, ph_dir;
			double inv_pdf, dir_pdf;
			c_light->SampleLightDir(rd, P_lit, N_lit, ph_dir, inv_pdf, dir_pdf);
			auto flux = N_lit.dot(ph_dir)
						* inv_pdf / (light_pmf * dir_pdf)
						* c_light->GetSurfaceProperty()->emission;

			assert(!(-1*flux).nonzero());

			// Cast Photon
			int castFlag = PM_DIRECT_LIGHT;
			if (p_caustic_ > -eps)
			{
				if (urand(rd) < 1 - p_caustic_)
				{
					flux *= 1 / (1 - p_caustic_);
					castFlag |= PM_NON_CAUSTIC;
				}
				else
				{
					flux *= 1 / p_caustic_;
					castFlag |= PM_CAUSTIC;
				}
			}
			else
				castFlag |= PM_CAUSTIC | PM_NON_CAUSTIC;

			CastPhoton(rd, P_lit, PhotonInfo{ph_dir, flux, 0}, pool, castFlag);

			// Update status & Merge photon samples
			if (++cnt % merge_every == 0)
			{
				_UpdatePMImage(pool);
				photons_.Merge(pool);
				pool.clear();
				cerr << "\rPhoton phase: " << photons_.vec.size() << " effective photons";
			}
		}
		_UpdatePMImage(pool);
		photons_.Merge(pool);
		pool.clear();
	};

	// Distribute threads
	for (int th = 0; th < n_threads; ++th)
	{
		int c_n_photons = n_photons / n_threads;
		if (th + 1 == n_threads)
			c_n_photons += n_photons % n_threads;

		threads_[th] = thread([this, &emitPhotons, c_n_photons, th]() {
			emitPhotons(rand_datum + th, c_n_photons);
			cerr << "Thread " << th << " finished\n";
		});
	}
	for (int th = 0; th < n_threads; ++th)
		threads_[th].join();

	clog << "\nPhotons emitted.\n";
}

void PM::_PreSampleRays (rand_data *rd, int n_threads, int n_rays)
{
	const double SPECULAR_RATIO = 0.7;
	int n_pixels = width_ * height_;
	int rays_per_pixel = (n_rays + n_pixels - 1) / n_pixels;

#ifndef SPECTRAL
	n_ray_sample_base_ = rays_per_pixel;
	n_ray_samples_ = new_matrix<int>(width_, height_, rays_per_pixel);
#else
	// Wavelength-dependent refraction needs far more ray samples.
	// We compensate for it here and redistribute samples.

	n_ray_samples_ = new_matrix<int>(width_, height_, 0);
	auto camera = scene_->GetCamera();

#pragma omp parallel for
	for (int x = 0; x < width_; ++x)
	{
		int th_id = omp_get_thread_num();
		auto rd = &rand_datum[th_id];
		for (int y = 0; y < height_; ++y)
		{
			for (int k = 0; k < rays_per_pixel || k < 10; ++k)
			{
				Ray ray;
				double dx = urand(rd), dy = urand(rd);
				camera->Project(rd, x + dx, y + dy, ray);

				ray.S += eps_mv * ray.D;
				Maybe<Collision> tmp = scene_->GetCollision(ray);

				if (tmp)
				{
					auto t = tmp.v.obj->GetSurfaceProperty()->SampleBxDF(rd);
					auto bxdf = dynamic_cast<SpectralSpecular *>(t->bxdf);
					n_ray_samples_[x][y] += int(bxdf != nullptr);
				}
			}
		}
		if (urand(rd) < 1e-2)
			clog << "\r Presampling rays: " << 100 * x / width_ << "% (" << th_id << ")";
	}
	clog << "\n";

	double sum = 1;
	for (int x = 0; x < width_; ++x)
		for (int y = 0; y < height_; ++y)
			sum += n_ray_samples_[x][y];

	n_ray_sample_base_ = int((n_rays * (1 - SPECULAR_RATIO) + n_pixels - 1) / n_pixels);
	for (int x = 0; x < width_; ++x)
		for (int y = 0; y < height_; ++y)
		{
			int delta = int(ceil(n_rays * SPECULAR_RATIO * n_ray_samples_[x][y] / sum));
			n_ray_samples_[x][y] = n_ray_sample_base_ + delta;
		}

#endif
}

int primes[61] = {
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
		83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
		191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283
};
inline int rev (const int i,const int p)
{
	if (i == 0) return i;
	return p - i;
}
double halton (const int b, int j)
{
	const int p = primes[b];
	double h = 0.0, f = 1.0 / (double) p, fct = f;
	while (j > 0)
	{
		h += rev(j % p, p) * fct;
		j /= p;
		fct *= f;
	}
	return h;
}

void PM::EmitRays (int n_threads)
{
	int w = width_, h = height_;

	vector<pair<int, int>> pblkque, pixque;
	for (int d = 0; d < w; d += 16)
	{
		for (int e = 0; e < h; e += 16)
			pblkque.push_back(make_pair(d, e));
	}
	GP::random_shuffle(pblkque.begin(), pblkque.end());
	for (auto p: pblkque)
	{
		for (int x = p.first; x < p.first + 16 && x < width_; ++x)
			for (int y = p.second; y < p.second + 16 && y < height_; ++y)
				pixque.push_back(make_pair(x, y));
	}

	using iter = decltype(pixque.begin());
	static int rs_counter{0};

	auto emitRays = [this](rand_data *rd, iter ps, iter pe) {
		auto camera = scene_->GetCamera();
		vector<PMRay> samples;
		for (auto pi = ps; pi != pe; ++pi)
		{
			int x = pi->first, y = pi->second;
			int n_samples = n_ray_samples_[x][y];
#ifdef SPECTRAL
			int halton_s = (rs_counter += n_samples) - n_samples + 30;
#endif
			for (int k = 0; k < n_samples; ++k)
			{
				Ray ray;
				double dx = urand(rd), dy = urand(rd);
				double lpdf = 1.;
				camera->Project(rd, x + dx, y + dy, ray);
				PMRay ret(&pixels_[x][y]);
				ret.sqrad = sqr(R0_);
				ret.weight = Spectrum(lpdf);
#ifdef SPECTRAL
				if (n_samples > n_ray_sample_base_)
				{
					int df = int(halton(0, halton_s + k) * (Spectrum::N)) % Spectrum::N;
					for (int d = 0; d < SampledSpectrum::N; ++d)
						ret.weight.vs[d] *= (d == df ? SampledSpectrum::N : 0);
				}
#endif
				if (CastRay(rd, ray, ret))
					samples.push_back(ret);
			}
			if (urand(rd) < 1e-5)
				cerr << "\r Effective rays: " << samples.size();
		}
		cerr << "\n";
		rays_.Merge(samples);
	};

	for (size_t th = 0, bs = pixque.size() / n_threads; th < n_threads; ++th)
	{
		auto s = pixque.begin() + bs * th;
		auto e = th + 1 == n_threads ? pixque.end() - 1 : s + bs;
		threads_[th] = std::thread([emitRays, this, th, s, e] () { emitRays(&rand_datum[th], s, e); });
	}
	for (int th = 0; th < n_threads; ++th)
		threads_[th].join();

	GP::random_shuffle(rays_.vec.begin(), rays_.vec.end());
}

bool PM::CastRay (rand_data *rd, const Ray &ray_, PMRay &res)
{
	Ray ray(ray_);
	for (int d = 0; d < ray_max_depth_; ++d)
	{
		ray.S += eps_mv * ray.D;
		Maybe<Collision> tmp = scene_->GetCollision(ray);
		if (!tmp)
			break;

		auto collision = tmp.v;
		const auto &surface = collision.obj->GetSurfaceProperty();
		auto bxdf_i = surface->SampleBxDF(rd);
		auto bxdf = bxdf_i->bxdf;
		Spectrum K = bxdf_i->GetColor(collision.obj, collision.pos);
		if (bxdf && bxdf->IsDirac())
		{
			Vec3 nD;
			double valid;
#ifndef SPECTRAL
			valid = bxdf->Sample(rd, ray.D, collision.norm, collision.inside, nD);
#else
			valid = dynamic_cast<SpectralSpecular*>(bxdf)
					->Sample(rd, ray.D, collision.norm, collision.inside, res.weight, nD);
#endif
			ray = Ray(collision.pos, nD);
			res.weight *= valid * K;
		}
		else
		{
			// We have chosen non-Dirac surface, or light
			// In both cases photon map should give correct result as long as
			// *direct light photons are stored*
			res.P = collision.pos;
			res.N = collision.norm;
			res.D = ray.D;
			res.bxdf = bxdf;
			if (bxdf != nullptr) // Use photon map to compute Ke
				res.weight *= K;

			assert(!((-1 * res.weight).nonzero()));
			return true;
		}
	}
	return false;
}

void PM::CastPhoton (rand_data *rd, Vec3 src, PhotonInfo photon, vector<pair<Vec3, PhotonInfo>> &pool,
					 int castMode)
{
#ifdef DIRECT_ONLY
	castMode = PM_DIRECT | PM_NON_CAUSTIC;
	d = 1;
#endif
	if (castMode & PM_NON_CAUSTIC)
		pool.push_back(make_pair(src, photon));

	int dmax = max_depth_;
	enum { LIGHT, DIFFUSE, SPECULAR }
	last_state = LIGHT;
	for (int d = 0; d < dmax; ++d)
	{
		src += eps_mv * photon.dir;
		Maybe<Collision> tmp = scene_->GetCollision(Ray(src, photon.dir));
		if (!tmp)
			break;

		auto collision = tmp.v;
		auto surface = collision.obj->GetSurfaceProperty();
		auto bxdf_i = surface->SampleBxDF(rd);
		auto bxdf = bxdf_i->bxdf;
		Spectrum K = bxdf_i->GetColor(collision.obj, collision.pos);
		PhotonInfo n_photon{Vec3{},
							photon.flux,
							photon.dis + collision.t};
		if (bxdf && bxdf->IsDirac())
		{
			double valid;
#ifndef SPECTRAL
			valid = bxdf->Sample(rd, photon.dir, collision.norm, collision.inside, n_photon.dir);
#else
			valid = dynamic_cast<SpectralSpecular*>(bxdf)
					->Sample(rd, photon.dir, collision.norm, collision.inside, n_photon.flux, n_photon.dir);
#endif
			n_photon.flux *= valid;
			n_photon.flux *= K;
			last_state = SPECULAR;
		}
		else
		{
			if (bxdf)
			{
				// Store photon
				bool store = false;
				store |= bool(castMode & PM_DIRECT_LIGHT) && (d == 0);
				store |= bool(castMode & PM_NON_CAUSTIC) && (last_state != SPECULAR);
				store |= bool(castMode & PM_CAUSTIC) && (last_state == SPECULAR);
				if (store)
					pool.push_back(make_pair(collision.pos, photon));

				double dir_pdf;
				dir_pdf = bxdf->Sample(rd, photon.dir, collision.norm, 0, n_photon.dir);
				Spectrum fw = K;
				auto k_fw = 1.0 / dir_pdf
					  * n_photon.dir.dot(collision.norm)
					  * bxdf->PDF(photon.dir, collision.norm, n_photon.dir);
				fw *= k_fw;
				n_photon.flux *= fw;
				last_state = DIFFUSE;
			}
			else // emission
				n_photon.flux = Spectrum(0.);
		}
		if (n_photon.flux.nonzero())
		{
			double p_roulette;
			if (d < min(dmax, 4))
				p_roulette = 1;
			else
				p_roulette = min(p_roulette_, n_photon.flux.y());

			if (urand(rd) < p_roulette)
			{
				n_photon.flux = 1.0 / p_roulette * n_photon.flux;
				src = collision.pos;
				photon = n_photon;
				continue;
			}
		}
		break;
	}
}

struct PPPMRayUpdater
{
	PMRay &r;
	PPPMRayUpdater (PMRay &ray): r(ray) {}
	void operator() (const Vec3 &pos, const PhotonInfo &info)
	{
		if (r.bxdf && info.dir.dot(r.N) > 0)
			// R0 too large. May confuse further assertions.
			return;

		double fK = r.bxdf == nullptr ? 1 :
					r.bxdf->PDF(r.D, r.N, -info.dir);
		r.tau += fK * info.flux;
		r.cnt_photons += 1;
	}
};

void PM::PPPM ()
{
	threads_ = vector<thread>((size_t)n_threads_);
	double R = R0_;

	_PreSampleRays(rand_datum, n_threads_, n_rays_);
	for (int round = 0; round < start_round_; ++round)
		R *= sqrt((round + PM_ALPHA) / (round + 1));

	for (int round = start_round_; round < n_rounds_; ++round)
	{
		total_photons_ = 0;
		rays_.Clear();
		EmitRays(n_threads_);
		photons_.Clear();
		EmitPhotons(n_threads_, photons_per_round_);
		PhotonTree<PhotonInfo> tree(&photons_.vec);

		Averager<double> kdt_qcnt(0);
		int counter = 0;
		typedef decltype(rays_.vec.begin()) iter;

		auto PPPMUpdateWorker = [&](iter s, iter e) {
			auto env_light = scene_->GetEnvLight();
			for (iter i = s; i < e; ++i) {
				++counter;
				if (counter % 10000 == 1)
				{
					cerr << "\rTree querying: " << counter << " ";
				}
				PMRay &ray = *i;
				PMPixel *pixel = ray.pixel;
				PPPMRayUpdater updateRay(ray);
				tree.Query(ray.P, R, updateRay);

				Spectrum radiance = ray.weight * ray.tau / (total_photons_ * pi * sqr(R));
				radiance += ray.weight * env_light;

				pixel->color.Add(radiance);
				image_->Set(pixel->scr_x, pixel->scr_y, pixel->color.Get());

				kdt_qcnt.Add(tree.__counter);
			}
		};

		size_t rays_per_thread = rays_.vec.size() / n_threads_;
		for (int i = 0; i < n_threads_; ++i)
		{
			auto s = rays_.vec.begin() + rays_per_thread * i;
			auto e = rays_.vec.begin() + rays_per_thread * (i + 1);
			threads_[i] = thread([s, e, &PPPMUpdateWorker](){ PPPMUpdateWorker(s, e); });
		}
		for (int i = 0; i < n_threads_; ++i)
			threads_[i].join();

		clog << "\nPhoton round " << round << " finished.\n";
		clog << "\t  Radius = " << R
		     << "\tMean KDTree query = " << kdt_qcnt.Get() << endl;
		clog << std::ios::fixed << std::setprecision(4)
		     << "C(" << _probe_x << ", " << _probe_y << ") = "
		     << pixels_[_probe_x][_probe_y].color.Get() << endl;
		image_->Dump(round);

		R *= sqrt((round + PM_ALPHA) / (round + 1));
	}
}


struct PPMRayUpdater
{
	PMRay &r;
	double PM_ALPHA;
	PPMRayUpdater (double alpha, PMRay &ray): PM_ALPHA(alpha), r(ray) {}
	void operator() (const Vec3 &pos, const PhotonInfo &info) {
		if ((pos - r.P).sqnorm() > r.sqrad || (r.bxdf && info.dir.dot(r.N) > 0))
			return;

		double fK = r.bxdf == nullptr ? 1 :
					r.bxdf->PDF(r.D, r.N, -info.dir);
		auto phi = fK * info.flux;
		double n_cnt_photons = 1.0 / PM_ALPHA + r.cnt_photons;
		r.cnt_photons += PM_ALPHA;
		double f = r.cnt_photons / n_cnt_photons;
		r.tau = (r.tau + phi) * f;
		r.sqrad *= f;
	}
};

void PM::PPM ()
{
	threads_ = vector<thread>((size_t)n_threads_);

	rays_.Clear();
	_PreSampleRays(&rand_datum[0], n_threads_, n_rays_);
	EmitRays(n_threads_);

	//for (int round = 0; round < n_rounds_; ++round)
	for (int round = start_round_; round < n_rounds_; ++round)
	{
		for (int x = 0; x < width_; ++x)
			for (int y = 0; y < height_; ++y)
				pixels_[x][y].color.Reset(Spectrum(0));

		photons_.Clear();
		EmitPhotons(n_threads_, photons_per_round_);

		PhotonTree<PhotonInfo> tree(&photons_.vec);

		Averager<double> ray_radius(0), kdt_qcnt(0);
		int counter = 0;
		typedef decltype(rays_.vec.begin()) iter;
		auto PPMUpdate = [&](rand_data *rd, iter s, iter e) {
			for (iter i = s; i < e; ++i) {
				++counter;
				if (urand(rd) < 0.001)
				{
					cerr << "\rTree querying: " << counter << " ";
				}
				PMRay &ray = *i;
				PMPixel *pixel = ray.pixel;
				PPMRayUpdater updateRay(PM_ALPHA, ray);
				tree.Query(ray.P, sqrt(ray.sqrad), updateRay);

				// Modifying initial shared radius does not affect correctness
				if (ray.cnt_photons < eps)
					ray.sqrad *= 1.0 / PM_ALPHA;

				Spectrum radiance = ray.weight * ray.tau / (total_photons_ * pi * ray.sqrad);
				radiance += ray.weight * scene_->GetEnvLight();

				pixel->color.Add(radiance);
				image_->Set(pixel->scr_x, pixel->scr_y, pixel->color.Get());

				kdt_qcnt.Add(tree.__counter);
				ray_radius.Add(sqrt(ray.sqrad));
			}
		};

		size_t rays_per_thread = rays_.vec.size() / n_threads_;
		for (int i = 0; i < n_threads_; ++i)
		{
			auto s = rays_.vec.begin() + rays_per_thread * i;
			auto e = rays_.vec.begin() + rays_per_thread * (i + 1);
			threads_[i] = thread([i, s, e, this, &PPMUpdate](){ PPMUpdate(&rand_datum[i], s, e); });
		}
		for (int i = 0; i < n_threads_; ++i)
			threads_[i].join();

		clog << "\nPhoton round " << round << " finished.\n";
		clog << "\t  Mean radius = " << ray_radius.Get()
		     << "\tMean KDTree query = " << kdt_qcnt.Get() << endl;
		clog << std::ios::fixed << std::setprecision(4)
		     << "C(" << _probe_x << ", " << _probe_y << ") = "
		     << pixels_[_probe_x][_probe_y].color.Get() << endl;
		image_->Dump(round);
	}

}

void PMParams::LoadFrom (const json &j)
{
#define L(d, X, v) X=J_AT_##d##_NZ(j, #X, v)
	L(I, max_depth, 11);
	L(I, ray_max_depth, 0);
	L(I, start_round, 0);
	L(F, R0, 0.02);
	L(F, p_caustic, 0.75);
	L(F, p_roulette, 0.98);
	L(F, alpha, 0.7);
	L(I, n_threads, 4);
	L(I, n_rays, 1e6);
	L(I, photons_per_round, 2e6);
	L(I, n_rounds, 100);
#ifndef NO_CV
	L(I, display, 1);
#else
	L(I, display, 0);
#endif
}
