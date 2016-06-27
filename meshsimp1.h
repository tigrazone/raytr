
//
// Created by dc on 5/27/16.
//

#ifndef RAYTR_MESHREDUCER1_H
#define RAYTR_MESHREDUCER1_H

#include <stdlib.h>
#include "common.h"
#include "image.h"

typedef Eigen::Matrix4d FMat;
typedef Eigen::Vector4d FVec;

const double epss = 1e-4;

inline Vec3 C (const FVec &u) { return Vec3(u(0), u(1), u(2)); }
inline FVec C (const Vec3 &u, double l)
{
	FVec ret;
	ret << u.x, u.y, u.z, l;
	return ret;
}

struct Face;

struct Vertex
{
	std::set<Face*> faces;
	std::set<Vertex*> vertices;
	FVec P;
	FMat V;
	bool removed, boundary;
	int hash;
	Vertex () { hash = 0; removed = boundary = false; }
	inline void add_face (Face *f);
	inline void replace (Vertex *u, Vertex *v);
	inline void update_V (int ctime);
};

struct Face
{
	Vertex *a, *b, *c;
	FMat m;
	Vec3 N_;
	int m_set;
	Face ()
	{
		a = b = c = nullptr;
		m_set = 0;
	}
	Face (Vertex *a_, Vertex *b_, Vertex *c_, const FMat &m_) :
			a(a_), b(b_), c(c_), m(m_), m_set(0) { }

	Vertex *operator[] (int id) { return (&a)[id]; }
	bool has (const Vertex *u) { return a == u || b == u || c == u; }
	void update_m ()
	{
#ifdef NO_UPDATE_M
		if (m_set)
			return;

		m_set = 1;
#endif
		assert (!(a->removed || b->removed || c->removed));
		Vec3 va = C(a->P);
		Vec3 vb = C(b->P);
		Vec3 vc = C(c->P);
		Vec3 N = (vb - va).normalized().cross(vc - va);
		if (N.sqnorm() > sqr(eps)) N_ = N = N.normalized();
		else N = N_;
		FVec t = C(N, -N.dot(va));
		m = t * t.transpose();
	}
	void replace (Vertex *u, Vertex *v)
	{
		assert((a==u)+(b==u)+(c==u)==1);
		assert((a==v)+(b==v)+(c==v)==0);
		for (int d = 0; d < 3; ++d) if ((&a)[d] == u) (&a)[d] = v;
		assert(!(a->P.isApprox(b->P, 1e-6)));
		assert(!(c->P.isApprox(b->P, 1e-6)));
		assert(!(a->P.isApprox(c->P, 1e-6)));
	}
};

void Vertex::add_face (Face *f)
{
	faces.insert(f);
	for (int k = 0; k < 3; ++k)
		if ((*f)[k] != this) vertices.insert((*f)[k]);
}

void Vertex::replace (Vertex *u, Vertex *v)
{
	assert(this != v);
	vertices.erase(u);
	vertices.insert(v);
}

class MeshSimplifier1
{
	vector<Vertex> vertices;
	vector<Face> faces;

public:
	MeshSimplifier1 (const string &objPath);
	void Run (double target_ratio);
};

void Vertex::update_V (int ctime)
{
	if (ctime == hash) return;
	hash = ctime;
	V *= 0;
	for (auto f: faces) V += f->m;
	if ((boundary = (vertices.size() > faces.size())))
	{
		FMat m;
		m << 1, 0, 0, -V(0),
			 0, 1, 0, -V(1),
			 0, 0, 1, -V(2),
			 -V(0), -V(1), -V(2), V(0) * V(0) + V(1) * V(1) + V(2) * V(2);
		V += 1e4 * m;
	}
}

#endif
