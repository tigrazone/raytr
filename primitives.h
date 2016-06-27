//
// Created by dc on 4/7/16.
//

#ifndef RAYTR_BASIC_OBJECTS_H
#define RAYTR_BASIC_OBJECTS_H

#include "objekt.h"
// FIXME: refactor this!
#include "mesh.h"

class Sphere: public Primitive
{
	Vec3 o_;
	double r_;
public:
	Sphere(Vec3 o, double r, PSurfaceProperty prop): Primitive(prop), o_(o), r_(r)
	{
		area_ = 4. / 3. * pi * r * r * r;
		bbox_ = BBox(o_ - Vec3(r, r, r), o_ + Vec3(r, r, r));
	}
	virtual Maybe<Collision> GetCollision (const Vec3 &src, const Vec3 &dir) const override;
	virtual void SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const override;
	virtual Vec3 SampleSurface (rand_data *rd, const Vec3 &src) const override;
	virtual double SAInvPDF (const Vec3 &src) const override;
	static Sphere* Create (const json &json_);

	virtual Vec3 ToUV (const Vec3 &surfacePoint) const override
	{ throw std::domain_error("not implemented"); }
};

class Triangle: public Primitive
{
protected:
	Vec3 A_, B_, C_;
	Vec3 N_, BmA_, CmB_, AmC_;
	CPVec3 NA_, NB_, NC_; // vertex normal. may be null (no normal interpolation)
	CPVec3 TA_, TB_, TC_; // texture coord. may be null
	double a_;

	// N_ points to the ``inside'' of triangle if it is part of a mesh

	void Init ()
	{
		BmA_ = B_ - A_;
		CmB_ = C_ - B_;
		AmC_ = A_ - C_;
		N_ = BmA_.normalized().cross(-AmC_).normalized(); // Mesh is small; avoid FP error
		a_ = N_.dot(A_);
		area_ = 0.5 * BmA_.cross(AmC_).norm();

		bbox_ = BBox::Min();
		for (int d = 0; d < 3; ++d)
		{
			for (int u = 0; u < 3; ++u)
			{
				bbox_.crd_min(d) = std::min(bbox_.crd_min(d), (&A_)[u][d]);
				bbox_.crd_max(d) = std::max(bbox_.crd_max(d), (&A_)[u][d]);
			}
		}
	}

public:
	Triangle (Vec3 A, Vec3 B, Vec3 C, PSurfaceProperty prop) :
			Primitive(prop), A_(A), B_(B), C_(C)
	{
		NA_ = NB_ = NC_ = TA_ = TB_ = TC_ = nullptr;
		Init();
	}

	Triangle (Vec3 A, Vec3 B, Vec3 C,
			  CPVec3 NA, CPVec3 NB, CPVec3 NC,
			  CPVec3 TA, CPVec3 TB, CPVec3 TC,
			  PSurfaceProperty prop) :
			Primitive(prop), A_(A), B_(B), C_(C),
			NA_(NA), NB_(NB), NC_(NC), TA_(TA), TB_(TB), TC_(TC)
	{
		Init();
		if (NA)
		{
			if (!((*NA).dot(N_) > -eps && (*NB).dot(N_) > -eps && (*NC).dot(N_) > -eps))
			{
				// may be hacking the model
				cerr << "Warning: incorrect surface normal loaded\n";
				NA_ = NB_ = NC_ = nullptr;
			}
		}
	}

	virtual Maybe<Collision> GetCollision (const Vec3 &src, const Vec3 &dir) const override;
	virtual void SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const override;
	virtual Vec3 SampleSurface (rand_data *rd, const Vec3 &src) const override;
	virtual double SAInvPDF (const Vec3 &src) const override;
	virtual Vec3 ToUV (const Vec3 &surfacePoint) const override;
	static Triangle *Create (const json &json_);
};

class Plane: public Primitive
{
protected:
	Vec3 N_;
	double a_; // x.dot(N) + a == 0; N.norm() == 1
public:
	Plane (Vec3 N, double a, PSurfaceProperty prop): Primitive(prop), N_(N), a_(a)
	{
		// TODO: ensure bbox is correctly handled.
		area_ = 0;
		bbox_ = BBox::Min();
	}
	virtual Maybe<Collision> GetCollision (const Vec3 &src, const Vec3 &dir) const override;
	virtual void SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const override;
	virtual Vec3 SampleSurface (rand_data *rd, const Vec3 &src) const override;
	virtual double SAInvPDF (const Vec3 &src) const override;
	static Plane* Create (const json &json_);
	Vec3 N () const { return N_; }

	virtual Vec3 ToUV (const Vec3 &surfacePoint) const override
	{ throw std::domain_error("not implemented"); }
};

class Rectangle: public Plane
{
	Vec3 A_, B_;
	double rA_[2], rB_[2];
	Vec3 P_[4]; // For solid angle calculation
public:
	Rectangle (const Plane &p,
			   Vec3 A, double rA1, double rA2,
			   Vec3 B, double rB1, double rB2):
			Plane(p), A_(A.normalized()), B_(B.normalized())
	{
		assert(rA1 <= rA2 && rB1 <= rB2);
		rA_[0] = rA1, rA_[1] = rA2;
		rB_[0] = rB1, rB_[1] = rB2;
		P_[0] = -a_ * N_ + rA2 * A_ + rB1 * B_;
		P_[1] = -a_ * N_ + rA1 * A_ + rB1 * B_;
		P_[2] = -a_ * N_ + rA2 * A_ + rB2 * B_;
		P_[3] = -a_ * N_ + rA1 * A_ + rB2 * B_;
		area_ = (rA2 - rA1) * (rB2 - rB1);
		bbox_ = BBox::Min();
		for (int d = 0; d < 3; ++d)
		{
			for (int u = 0; u < 4; ++u)
			{
				bbox_.crd_min(d) = std::min(bbox_.crd_min(d), P_[u][d]);
				bbox_.crd_max(d) = std::max(bbox_.crd_max(d), P_[u][d]);
			}
		}
	}
	virtual Maybe<Collision> GetCollision (const Vec3 &src, const Vec3 &dir) const override;
	virtual void SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const override;
	virtual Vec3 SampleSurface (rand_data *rd, const Vec3 &src) const override;
	virtual double SAInvPDF (const Vec3 &src) const override;
	static Rectangle* Create (const json &json_);

	virtual Vec3 ToUV (const Vec3 &surfacePoint) const override
	{ throw std::domain_error("not implemented"); }
};

inline Primitive* CreateBasicObject (const json &json_)
{
	std::string type = json_.at("t").get<std::string>();
	Primitive *ret;
	if (type == "sphere")
	{
		ret = Sphere::Create(json_);
	}
	else if (type == "triangle" || type == "tri")
	{
		ret = Triangle::Create(json_);
	}
	else if (type == "plane")
	{
		// Planes are deprecated.
		Plane *p = Plane::Create(json_);
		Vec3 A, B;
		GetPerpVecs(p->N(), A, B);
		ret = new Rectangle(*p, A, -1e3, 1e3, B, -1e3, 1e3);
	}
	else if (type == "rectangle" || type == "rect")
	{
		ret = Rectangle::Create(json_);
	}
	else if (type == "mesh")
	{
		ret = Mesh::Create(json_);
	}
	else
		assert(false);

	if (json_.count("ldrange"))
	{
		ret->SetUnidirectedLight(
				Vec3::Parse(json_["ldrange"][0]),
				Vec3::Parse(json_["ldrange"][1]));
	}
	return ret;
}

#endif //RAYTR_BASIC_OBJECTS_H
