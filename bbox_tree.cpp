//
// Created by dc on 5/7/16.
//

#include "bbox_tree.h"

void BBoxTree::Query_ (uint32_t node, const Vec3 &S, const Vec3 &D,
					   double t0, double t1,
					   double &res_t, const Primitive *&res_obj) const
{
	if ((B[node].flag & 3) == 3) // LEAF
	{
		uint32_t s_prim = B[node].p_prims;
		uint32_t e_prim = s_prim + (B[node].n_prims >> 2u);
		res_t = inf;
		for (auto j = s_prim; j < e_prim; ++j)
		{
			auto r = prims[j]->GetCollision(S, D);
			if (r.valid && r.v.t < res_t)
			{
				res_t = r.v.t;
				res_obj = r.v.obj;
			}
		}
	}
	else
	{
		int ax = B[node].flag & 3;
		double inv_d = fabs(D[ax]) < eps ? 1.0 / eps : 1.0 / D[ax];
		double t = (B[node].split - S[ax]) * inv_d;
		// if D[ax] ~ 0, split > S[ax] ? vis lch : vis rch

		// ch[0]: nearer; ch[1]: farther
		uint32_t ch[2] = {node + 1u, B[node].p_rch >> 2u};
		if (inv_d < 0) swap(ch[0], ch[1]);

		if (t - eps > t1)
			Query_(ch[0], S, D, t0, t1, res_t, res_obj);
		else if (t + eps < t0)
			Query_(ch[1], S, D, t0, t1, res_t, res_obj);
		else
		{
			double rt;
			const Primitive *robj;
			Query_(ch[0], S, D, t0, t, rt, robj);
			if (rt > t)
				Query_(ch[1], S, D, t, t1, rt, robj);

			res_t = rt;
			res_obj = robj;
		}
	}
}

struct SAHEvent
{
	static const int ENTER = -1, LEAVE = 1;
	double x;
	int t, id;
	bool operator< (const SAHEvent &b) const
	{
		return fabs(x - b.x) < eps ? t < b.t : x < b.x;
	}
};

static int SAH_SetupEvents (PCPrimitive const *prims, int n_prims, int axis, vector<SAHEvent> &events)
{
	int n_events = 2 * n_prims;
	events.resize(n_events);
#pragma omp parallel for
	for (int i = 0; i < n_prims; ++i)
	{
		const auto &bbox = prims[i]->GetBBox();
		events[i * 2] = {bbox.crd_min(axis), SAHEvent::ENTER, i};
		events[i * 2 + 1] = {bbox.crd_max(axis), SAHEvent::LEAVE, i};
	}
	GP::sort(events.begin(), events.end());
	return n_events;
}

static double SAH_Split (PCPrimitive const *prims, int n_prims, int axis, BBox gbox, float &split)
{
	// Create & sort scanning line events
	vector<SAHEvent> events;
	int n_events;
	n_events = SAH_SetupEvents(prims, n_prims, axis, events);

	BBox rbox = gbox, lbox = gbox;
	double best_cost = inf;
	for (int j = 0, e, n_overlap = 0; j < n_events; )
	{
		// Consider splitting at events[j].pos + deps. Split events into [0, i) and [i, n).
		for (e = j; e < n_events && fabs(events[e].x - events[j].x) < eps; ++e) ;

		// Calculate cost
		for (int k = j; k < e; ++k)
			n_overlap -= events[k].t;

		float c_split = (float)events[j].x;
		c_split += (e == n_events) ? 16 * eps : 1e-3 * (events[e].x - events[j].x);

		lbox.crd_max(axis) = c_split;
		rbox.crd_min(axis) = c_split;

		int n_l = (e - n_overlap) / 2 + n_overlap,
			n_r = (n_events - e - n_overlap) / 2 + n_overlap;
		double cost = lbox.SA() / gbox.SA() * n_l +
				rbox.SA() / gbox.SA() * n_r;

		if (cost < best_cost)
		{
			best_cost = cost;
			split = c_split;
		}

		j = e;
	}

	return best_cost;
}

uint32_t BBoxTree::Build_ (uint32_t node, PCPrimitive const *prims, uint32_t n_prims, BBox bbox, int depth)
{
	const double Ki = 80; // Cost(isect) / Cost(call)

	// Allocate storage for node
	assert(B_container_.size() == node);
	B_container_.push_back(BTNode{});

	// Find best splitting cost
	double spl_cost = inf;
	double leaf_cost = n_prims * Ki;
	float spl_pos;
	int spl_axis;
	int preferred_axis = bbox.max_axis();
	for (int d = 0; d < 3; ++d)
	{
		int ax = (preferred_axis + d) % 3;
		float c_sp;
		double c_co = SAH_Split(prims, n_prims, ax, bbox, c_sp) * Ki + 1;
		if (c_co < spl_cost)
		{
			spl_cost = c_co;
			spl_pos = c_sp;
			spl_axis = ax;
		}
		if (ax == preferred_axis && spl_cost < leaf_cost)
			break;
	}

	if (spl_cost < leaf_cost && depth < MAX_DEPTH)
	{
		// CREATE INTERIOR NODE
		vector<SAHEvent> events;
		vector<PCPrimitive> pbuff(n_prims);
		int n_events, cn_prims = 0;

		n_events = SAH_SetupEvents(prims, n_prims, spl_axis, events);

		// Select (unique) primitives which intersect left box
		for (int d = 0; d < n_events; ++d)
		{
			if (events[d].x < spl_pos + eps && events[d].t == SAHEvent::ENTER)
				pbuff[cn_prims++] = prims[events[d].id];
		}

		// Build left subtree
		uint32_t next_avail_pos = node + 1;
		BBox cbox(bbox);
		cbox.crd_max(spl_axis) = spl_pos;
		next_avail_pos = Build_(next_avail_pos, &pbuff[0], cn_prims, cbox, depth + 1);

		// Now we know everything to init node
		B_container_[node].InitInterior(next_avail_pos, spl_axis, spl_pos);

		// Select unique primitives for rch
		cn_prims = 0;
		for (int d = 0; d < n_events; ++d)
		{
			if (spl_pos < events[d].x + eps && events[d].t == SAHEvent::LEAVE)
				pbuff[cn_prims++] = prims[events[d].id];
		}

		// Build rch
		cbox = bbox;
		cbox.crd_min(spl_axis) = spl_pos;
		next_avail_pos = Build_(next_avail_pos, &pbuff[0], cn_prims, cbox, depth + 1);

		return next_avail_pos;
	}
	else
	{
		// INITIALIZE AS LEAF

		// Allocate permanent storage for primitives
		prims_container_.reserve(prims_container_.size() + n_prims);
		uint32_t p_prim = (uint32_t)prims_container_.size();
		for (int u = 0; u < n_prims; ++u)
			prims_container_.push_back(prims[u]);

		// Init
		B_container_[node].InitLeaf(n_prims, p_prim);
		return node + 1u;
	}
}