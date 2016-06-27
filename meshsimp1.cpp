//
// Created by dc on 5/30/16.
//

#include "meshsimp1.h"
#include "objloader.h"
#include <tiny_obj_loader.h>
#include <queue>
#include <fstream>
#include <unordered_map>
using namespace std;

typedef unsigned long long ULL;

MeshSimplifier1::MeshSimplifier1 (const string &objPath)
{
	using namespace tinyobj;

	vector<shape_t> shapes;
	vector<material_t> materials;

	string err;
	bool ret = LoadObj(shapes, materials, err, objPath.c_str());

	if (!err.empty())
	{
		// `err` may contain warning message.
		cerr << "When loading mesh: " << err << endl;
	}

	if (!ret)
		throw logic_error("Mesh failed to load");

	cout << "# of shapes    : " << shapes.size() << endl;
	cout << "# of materials : " << materials.size() << endl;

	for (size_t i = 0; i < shapes.size(); ++i)
	{
		const auto &m_indices = shapes[i].mesh.indices;
		this->vertices.insert(vertices.end(), shapes[i].mesh.positions.size() / 3,
							  Vertex());
		this->faces.insert(faces.end(), m_indices.size() / 3, Face());
	}

	auto vpos = vertices.begin();
	auto fpos = faces.begin();

	rand_data rd;
	rand_init(&rd);

	for (size_t i = 0; i < shapes.size(); ++i)
	{
		printf("shape[%ld].vertices: %ld\n", i, shapes[i].mesh.positions.size());
		vector<Vertex*> m_verts;
		for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++)
		{
			FVec vv;
			vv << shapes[i].mesh.positions.at(3 * v + 0),
					shapes[i].mesh.positions.at(3 * v + 1),
					shapes[i].mesh.positions.at(3 * v + 2),
					1;
			vpos->P = vv;
			m_verts.push_back(&*vpos);
			vpos++;
		}

		const auto &m_indices = shapes[i].mesh.indices;
		for (size_t f = 0; f < m_indices.size() / 3; f++)
		{
			int a = m_indices[3 * f + 2], b = m_indices[3 * f + 1], c = m_indices[3 * f + 0];
			// Face
			*fpos = Face{ m_verts[a], m_verts[b], m_verts[c], FMat() };
			fpos->update_m();

			for (int j = 0; j < 3; ++j)
				fpos[0][j]->add_face(&*fpos);

			fpos++;
		}
	}
}

typedef tuple<double, FVec, Vertex*, Vertex*, int> qnode;
inline bool QNode (Vertex *x, Vertex *y, int time, qnode &res)
{
	static atomic<int> counter{ 0 };
	x->update_V(time);
	y->update_V(time);

	if (x->boundary && y->boundary)
		return false;

	FMat m = x->V + y->V;
	m(3, 0) = m(3, 1) = m(3, 2) = 0; m(3, 3) = 1;
	FVec v; v << 0, 0, 0, 1;

//	FVec v1 = m.colPivHouseholderQr().solve(v);
	FVec v1 = m.householderQr().solve(v);
	// normal equation (ldlt) gives incorrect result

	if (!(m * v1).isApprox(v, 1e-3))
	{
		v1 = (x->P + y->P) / 2;
	}
	double cost = (v1.transpose() * (x->V + y->V) * v1);
	assert(!std::isnan(cost));
	res = qnode(cost, v1, x, y, time);
	return true;
}

void MeshSimplifier1::Run (double target_ratio)
{
#define IDV(p) (p-&vertices[0])
#define UPAIR(x, y) ((IDV(min(x,y)) << 32ULL) | IDV(max(x,y)))
	struct qcomp
	{
		bool operator() (const qnode &a, const qnode &b) { return get<0>(a) > get<0>(b); }
	};
	priority_queue<qnode, vector<qnode>, qcomp> pque;
	unordered_map<ULL, int> last_updated;

	int ie = vertices.size();
	for (int i = 0; i < ie; i++)
	{
		if (i % 100 == 0)
			cerr << "\r Iter" << i;

		auto u = &vertices[i];
		for (auto v: u->vertices)
			{
				qnode q;
				if (! QNode(u, v, 1, q))
					continue;

				pque.push(q);
				last_updated[UPAIR(u, v)] = 1;
			}
	} cerr << endl;

	size_t n_iter = faces.size() * (1 - target_ratio);
	double total_cost = 0;
	for (size_t rd = 1; rd <= n_iter; )
	{
		double cost; FVec up1; Vertex *u, *v; int times;
		tie(cost, up1, u, v, times) = pque.top(); pque.pop();
		if (last_updated[UPAIR(u, v)] != times)
			continue;

		if (u->removed || v->removed) continue;

#ifdef NDEBUG
		if (rd % 1000 == 0)
#endif
		{
			cerr << "Round " << rd << endl;
			cerr << u - &vertices[0] << ' ' << v - &vertices[0] << ' ' << cost << endl;
		}

#ifndef NDEBUG
		{
			for (auto &f: faces)
				if (!f.a->removed && !f.b->removed && !f.c->removed)
					f.update_m();

			qnode q;
			bool r = QNode(u, v, -1, q);
			auto exp = get<0>(q);
			assert(fabs((exp - cost) / cost) < 1e-3);
		}
#endif

		total_cost += cost;
		v->P = u->P = up1;
		v->removed = true;

		// Clean faces
		for (auto i = v->faces.begin(); i != v->faces.end(); )
		{
			auto f = *i;
			assert(f->has(v));
			if (!u->faces.count(f))
			{
				// Rename v <= u;
				f->replace(v, u);
				u->faces.insert(f);
				i++;
			}
			else
			{
				// Remove f
				++rd;
				auto ni = next(i);
				f->a->faces.erase(f);
				f->b->faces.erase(f);
				f->c->faces.erase(f);
				i = ni;
			}
		}

		// Vertices
		for (auto vv: v->vertices)
			if (!vv->removed && vv != u)
			{
				vv->replace(v, u);
				u->vertices.insert(vv);
			}

		u->vertices.erase(v);
		for (auto x: u->vertices)
			x->removed = x->removed || x->vertices.empty();

		for (auto f: u->faces)
			f->update_m();

		// Update (edge) costs
		for (auto x: u->vertices)
			for (auto y: x->vertices)
			{
				assert(x->vertices.count(u));
				if (!x->removed && !y->removed)
				{
					auto cx = x, cy = y;
					qnode q;
					if (!QNode(cx, cy, rd, q))
						continue;

					if (last_updated[UPAIR(cx, cy)] != rd)
					{
						pque.push(q);
						last_updated[UPAIR(cx, cy)] = rd;
					}
				}
			}
	}

	// Dump obj

	ofstream fout("dumped.obj");
	fout << "o Object\n";
	for (int d = 0; d < vertices.size(); ++d)
	{
		fout << "v " << vertices[d].P(0) << " " << vertices[d].P(1) << " " << vertices[d].P(2) << endl;
	}
	for (auto &f: faces)
	{
		if (!f.a->removed && !f.b->removed && !f.c->removed)
			fout << "f " << 1 + IDV(f.c) << ' ' << 1 + IDV(f.b) << ' ' << 1 + IDV(f.a) << endl;
	}
	fout.close();
}
