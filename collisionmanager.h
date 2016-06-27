//
// Created by dc on 5/7/16.
//

#ifndef RAYTR_COLLISIONMANAGER_H
#define RAYTR_COLLISIONMANAGER_H

#include "objekt.h"
// #include "mesh.h"

class CollisionManager
{
public:
	virtual Maybe<Collision> Get (const Ray &r) const = 0;
};

class BFCollisionMGR: public CollisionManager
{
	std::vector<const Primitive*> primitives_;
public:
	BFCollisionMGR (std::vector<Primitive*> prims)
	{
		for (auto p: prims) primitives_.push_back(p);
	}
	virtual Maybe<Collision> Get (const Ray &r) const override
	{
		auto ret = Maybe<Collision>::None();
		for (const auto &obj: primitives_)
		{
			Maybe<Collision> cur = obj->GetCollision(r.S, r.D);
			if (cur && (!ret || ret.v.t > cur.v.t))
				ret = cur;
		}
		return ret;
	}
};

#endif //RAYTR_COLLISIONMANAGER_H
