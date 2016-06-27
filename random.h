//
// Created by dc on 4/9/16.
//

#ifndef RAYTR_RANDOM_H
#define RAYTR_RANDOM_H

#include <random>
#include <dSFMT/dSFMT.h>
#include <cstdlib>
#include <ctime>

#include "common.h"

#define VALGRIND

class Random
{
	Random ()
	{
		srand(time(NULL));
	}
public:
	static Random random;
};

typedef dsfmt_t rand_data;

inline void rand_init (rand_data *r)
{
	dsfmt_init_gen_rand(r, rand());
}
inline double urand (rand_data *r)
{
	return dsfmt_genrand_close_open(r);
}
inline void urand_arr (rand_data *r, double *arr, int sz)
{
	dsfmt_fill_array_open_close(r, arr, sz);
}
inline double urand (rand_data *r, double lo, double hi) { return urand(r) * (hi - lo) + lo; }
inline double urand (rand_data *r, double hi) { return urand(r) * hi; }

inline Vec3 RandomRotate (rand_data *rd, double r, const Vec3 &N)
{
	double theta = urand(rd, 2 * pi);
	double x = r * cos(theta), y = r * sin(theta),
			z = sqrt(std::max(0.0, 1 - sqr(r)));
	// Rotate z |-> N
	Vec3 A, B;
	GetPerpVecs(N, A, B);
	Vec3 ret = x * A + y * B + z * N;
	return ret;
}

inline Vec3 CosRandHemisphere (rand_data *rd, const Vec3 &N)
{
	return RandomRotate(rd, sqrt(urand(rd)), N);
}
inline void CosRandHemisphere (rand_data *rd, const Vec3 &N, Vec3 &res, double &pdf)
{
	res = CosRandHemisphere(rd, N);
	pdf = res.dot(N) / pi;
}

inline void URandPlate (rand_data *rd, double r, double &x, double &y)
{
	do {
		x = urand(rd) * 2 - 1;
		y = urand(rd) * 2 - 1;
	} while (sqr(x) + sqr(y) > 1);
	x *= r;
	y *= r;
}

#endif //RAYTR_RANDOM_H
