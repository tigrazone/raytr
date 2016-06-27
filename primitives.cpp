//
// Created by dc on 4/8/16.
//

#include "primitives.h"
#include "shade.h"
#include "random.h"

Sphere* Sphere::Create (const json &j)
{
	return new Sphere(Vec3::Parse(j.at("O")),
					  J_AT_F(j, "r"),
					  ShadingProperty::Create(j.at("surface")));
}

Plane* Plane::Create (const json &j)
{
	return new Plane(Vec3::Parse(j.at("N")).normalized(),
					 J_AT_F(j, "a"),
					 ShadingProperty::Create(j.at("surface")));
}

Rectangle* Rectangle::Create (const json &j)
{
	Plane *plane = Plane::Create(j);
	auto A = j.at("A");
	auto B = j.at("B");
	Rectangle *rect = new Rectangle(
			*plane,
			Vec3::Parse(A[0]), A[1].get<double>(), A[2].get<double>(),
			Vec3::Parse(B[0]), B[1].get<double>(), B[2].get<double>());
	delete plane;
	return rect;
}

Triangle* Triangle::Create (const json &j)
{
	return new Triangle(
			Vec3::Parse(j.at("A")),
			Vec3::Parse(j.at("B")),
			Vec3::Parse(j.at("C")),
			ShadingProperty::Create(j.at("surface")));
}

Maybe<Collision> Sphere::GetCollision (const Vec3 &src, const Vec3 &dir) const
{
	Vec3 l = o_ - src;
	double t_p = l.dot(dir),
			sqdis = l.dot(l) - t_p * t_p;
	if (sqdis > sqr(r_) + eps)
		return Maybe<Collision>::None();

	double t_ = sqrt(std::max(0., sqr(r_) - sqdis)),
			t;

	if (l.sqnorm() + eps > sqr(r_) && t_p < eps)
		return Maybe<Collision>::None();

	bool toInside = true;
	if ((o_ - src).sqnorm() > sqr(r_))
		t = t_p - t_;
	else
	{
		t = t_p + t_;
		toInside = false;
	}
	Vec3 pos = src + t * dir,
			N_c = (toInside ? 1 : -1) * (pos - o_).normalized();
	return Maybe<Collision>::Some({t, pos, N_c, this, toInside});
}

Maybe<Collision> Triangle::GetCollision (const Vec3 &src, const Vec3 &dir) const
{
	double d = N_.dot(-dir); // defn of N is opposite to Plane. fix that.
	if (fabs(d) > eps)
	{
		double t = -(a_ - N_.dot(src)) / d;
		if (t > eps && t < 1e5)
		{
			auto X = src + t * dir;
			auto t1 = BmA_.cross(X - A_).dot(N_);
			auto t2 = CmB_.cross(X - B_).dot(N_);
			auto t3 = AmC_.cross(X - C_).dot(N_);
			if (t1 < -eps || t2 * 4 < eps || t3 * 4 < eps) // avoid duplication
				return Maybe<Collision>::None();

			auto cur = !NA_ ? N_ : (t1 * *NC_ + t2 * *NA_ + t3 * *NB_).normalized();
			auto norm = d > eps ? cur : -cur;
			return Maybe<Collision>::Some({t, X, norm, this, d <= eps});
		}
	}
	return Maybe<Collision>::None();
}

Maybe<Collision> Plane::GetCollision (const Vec3 &src, const Vec3 &dir) const
{
	double d = N_.dot(dir);
	if (fabs(d) > eps)
	{
		double t = -(a_ + N_.dot(src)) / d;
		if (t > eps && t < 1e5)
		{
			Vec3 dst = src + t * dir,
					norm = d > eps ? -N_ : N_;
			// For plane "inside" is defined as the opposite direction of N_.
			return Maybe<Collision>::Some({t, dst, norm, this, d <= eps});
		}
	}
	return Maybe<Collision>::None();
}

Maybe<Collision> Rectangle::GetCollision (const Vec3 &src, const Vec3 &dir) const
{
	Maybe<Collision> col_plane = Plane::GetCollision(src, dir);
	if (col_plane)
	{
		Vec3 v = col_plane.v.pos;
		double dA = v.dot(A_), dB = v.dot(B_);
		if (dA + eps < rA_[0] || rA_[1] + eps < dA ||
			dB + eps < rB_[0] || rB_[1] + eps < dB)
			return Maybe<Collision>::None();
	}
	return col_plane;
}

void Plane::SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const
{
	throw std::domain_error("Plane cannot be light source.");
}
Vec3 Plane::SampleSurface (rand_data *rd, const Vec3 &) const
{
	throw std::domain_error("Plane cannot be light source.");
}
double Plane::SAInvPDF (const Vec3 &) const
{
	throw std::domain_error("Plane cannot be light source.");
}

inline double SolidAngle (const Vec3 A[], const Vec3 &b)
{
	Vec3 R[3] = {A[0] - b, A[1] - b, A[2] - b};
	double R0 = R[0].norm(), R1 = R[1].norm(), R2 = R[2].norm();
	double num = R[0].dot(R[1].cross(R[2])),
		   den = R0 * R1 * R2 + R[0].dot(R[1]) * R2 + R[0].dot(R[2]) * R1 + R[1].dot(R[2]) * R0;
	return 2 * fabs(atan2(num, den));
}

void Triangle::SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const
{
	double r1, r2;
	do {
		r1 = urand(rd), r2 = urand(rd);
	} while (r1 + r2 > 1);
	sample = (1 - r1 - r2) * A_ + r1 * B_ + r2 * C_;
//	double sr1 = sqrt(r1);
//	double a = 1 - sr1, b = sr1 * (1 - r2), c = r2 * sr1;
//	sample = a * A_ + b * B_ + c * C_;
	N = (urand(rd) < 0.5 ? 1 : -1) * N_;
}

Vec3 Triangle::SampleSurface (rand_data *rd, const Vec3 &src) const
{
	double r1, r2;
	do {
		r1 = urand(rd), r2 = urand(rd);
	} while (r1 + r2 > 1);
	return (1 - r1 - r2) * A_ + r1 * B_ + r2 * C_;
}

double Triangle::SAInvPDF (const Vec3 &src) const
{
	return SolidAngle(&A_, src);
}

void Rectangle::SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const
{
	sample = -a_ * N_ + urand(rd, rA_[0], rA_[1]) * A_ + urand(rd, rB_[0], rB_[1]) * B_;
	N = (urand(rd) < 0.5 ? 1 : -1) * N_;
}

Vec3 Rectangle::SampleSurface (rand_data *rd, const Vec3 &src) const
{
	return -a_ * N_ + urand(rd, rA_[0], rA_[1]) * A_ + urand(rd, rB_[0], rB_[1]) * B_;
}

double Rectangle::SAInvPDF (const Vec3 &src) const
{
	return SolidAngle(P_, src) + SolidAngle(P_ + 1, src);
}

void Sphere::SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const
{
	double phi = urand(rd, 2 * pi), theta = urand(rd, pi);
	double costh = cos(theta), sinth = sqrt(std::max(0., 1 - sqr(costh)));
	N = Vec3{sinth * cos(phi), sinth * sin(phi), costh};
	sample = o_ + r_ * N;
}

Vec3 Sphere::SampleSurface (rand_data *rd, const Vec3 &src) const
{
	Vec3 A, B, N = (src - o_).normalized();
	GetPerpVecs(N, A, B);
	double sinva = r_ / (src - o_).norm(),
		   va = sinva + eps > 1 ? 0 : asin(r_ / (src - o_).norm());
	double phi = urand(rd, 2 * pi), theta = urand(rd, va + eps, pi - va - eps),
			costh = cos(theta),
			x = r_ * costh * cos(phi),
			y = r_ * costh * sin(phi),
			z = r_ * sqrt(std::max(0., 1 - sqr(costh)));
	Vec3 sample = o_ + x * A + y * B + z * N;
	return sample;
}

double Sphere::SAInvPDF (const Vec3 &src) const
{
	// InvPDF
	double sqSinVA = sqr(r_) / (src - o_).sqnorm(),
		   cosVA = sqrt(std::max(0., 1 - sqSinVA)),
		   invq = 2 * pi * (1.0 - cosVA);
	// When src is inside sphere this gives 1/2pi as expected.
	return invq;
}

Vec3 Triangle::ToUV (const Vec3 &X) const
{
	if (!TA_)
		throw std::domain_error("not supported");

	auto tC = BmA_.cross(X - A_).dot(N_);
	auto tA = CmB_.cross(X - B_).dot(N_);
	auto tB = AmC_.cross(X - C_).dot(N_);
	return (tC * *TC_ + tB * *TB_ + tA * *TA_) / (tA + tB + tC);
}
