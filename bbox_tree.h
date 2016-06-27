//
// Created by dc on 5/7/16.
//

#ifndef RAYTR_BBOX_TREE_H
#define RAYTR_BBOX_TREE_H

#include "objekt.h"
#include "collisionmanager.h"

struct BTNode
{
	static const int LEAF = 3;
	// Lower 2 bits: 3 - leaf; 0 ~ 2 - interior (indicating axis to split)
	// Upper 30 bits: #Primitives
	union
	{
		uint32_t flag;
		uint32_t n_prims;
		uint32_t p_rch;
	};
	union
	{
		float split;
		uint32_t p_prims;
	};
	inline void InitInterior (uint32_t p_rch, int axis, float split)
	{
		this->split = split;
		this->flag = (p_rch << 2u) | (uint32_t)axis;
	}
	inline void InitLeaf (uint32_t n_prims, uint32_t p_prims)
	{
		this->flag = (n_prims << 2u) | 3u;
		this->p_prims = p_prims;
	}
};

using PCPrimitive = const Primitive*;

class BBoxTree: public CollisionManager
{
private:
	vector<BTNode> B_container_;
	vector<PCPrimitive> prims_container_;
	vector<PCPrimitive> all_prims_;

protected:
	static const int MAX_DEPTH = 20;
	BTNode *B;
	PCPrimitive *prims;
	BBox glob_;

	void Query_ (uint32_t node, const Vec3 &S, const Vec3 &D,
				 double t0, double t1,
				 double &res_t, const Primitive *&res_obj) const;
	uint32_t Build_ (uint32_t node, PCPrimitive const *prims, uint32_t n_prims, BBox bbox, int depth);

public:
	BBoxTree (vector<Primitive*> primitives)
	{
		for (auto p: primitives) all_prims_.push_back(p);
		glob_ = BBox::Min();
		for (auto p: all_prims_) glob_.merge_with(p->GetBBox());
		Build_(0, &all_prims_[0], (uint32_t)all_prims_.size(), glob_, 0);
		B = &B_container_[0];
		prims = &prims_container_[0];
		// all_prims_.clear();
	}

	virtual Maybe<Collision> Get (const Ray &r) const override
	{
		double rt;
		const Primitive *robj;
		Maybe<Collision> ret;
		Query_(0, r.S, r.D, 0, inf, rt, robj);
		if (rt < inf * 0.5)
			ret = robj->GetCollision(r.S, r.D);
		else
			ret = Maybe<Collision>::None();

		return ret;
	}
};

class BBTChecker: public CollisionManager
{
private:
	BBoxTree *bt;
	BFCollisionMGR *mgr;
public:
	BBTChecker (vector<Primitive*> prims)
	{
		cerr << prims.size();
		bt = new BBoxTree(prims);
		mgr = new BFCollisionMGR(prims);
	}

	virtual Maybe<Collision> Get (const Ray &r) const override
	{
		auto btret = bt->Get(r);
		auto bfret = mgr->Get(r);
		if (btret.valid != bfret.valid ||
			btret.valid && btret.v.t - 1e-5 > bfret.v.t)
		{
			cerr << r.S << " -> " << r.D << ": " <<
			"act: " << btret.valid << " " << btret.v.t << " " << btret.v.pos << ";\t" <<
			"exp: " << bfret.valid << " " << bfret.v.t << " " << bfret.v.pos << endl;
			throw std::logic_error("Ooooooops");
		}
		return btret;
	}


};

#endif //RAYTR_BBOX_TREE_H
