//
// Created by dc on 4/7/16.
//

#ifndef RAYTR_SCENE_H
#define RAYTR_SCENE_H

#include "common.h"
#include "random.h"
#include "shade.h"
#include "camera.h"
#include "collisionmanager.h"
#include <string>
#include <list>

struct PointLight
{
	Vec3 O;
	Spectrum intensity;
	PointLight (Vec3 O_, Spectrum intensity_): O(O_), intensity(intensity_) {}
	static PointLight* Create (const json &json_);
};

class Scene
{
	vector<Primitive*> objects_; // Stores all objects
	// Stores all lights. \subseteq objects_; second = light pmf
	vector<pair<Primitive*, double>> lights_;
	Spectrum env_light_;

	vector<Destroyable*> resources_;
	CollisionManager *collision_;
	Camera *camera_;
	int width_;
	int height_;

public:
	Scene (vector<Primitive*> objects,
		   vector<Primitive*> lights,
		   vector<Destroyable*> resources,
		   CollisionManager *collision,
		   Camera *camera,
		   Spectrum env_light,
		   int width = 800,
		   int height = 600):
			objects_(objects), resources_(resources),
			collision_(collision),
			camera_(camera),
			env_light_(env_light),
			width_(width), height_(height)
	{
		// Calculate Light PMF
		lights_.resize(lights.size());
		std::transform(lights.begin(), lights.end(), lights_.begin(), [](Primitive * obj){
			auto c = obj->GetSurfaceProperty()->emission;
			double w = 0;
			for (int d = 0; d < Spectrum::N; ++d)
				w += sqr(c[d]);

			return make_pair(obj, sqrt(w));
		});
		double sum_power = 0;
		for (auto &i: lights_)
			sum_power += i.second;

		for (auto &i: lights_)
			i.second /= sum_power;
	}
	Scene ()
	{
		for (auto obj: objects_)
		{
			delete obj;
		}
		for (auto res: resources_)
		{
			delete res;
		}
		delete camera_;
		delete collision_;
	}

	const Camera* GetCamera () const { return camera_; }
	const vector<pair<Primitive*, double>> & GetLights () const { return lights_; }
	Spectrum GetEnvLight () const { return env_light_; }

	const Primitive *SampleLight (rand_data *rd, double &pmf) const
	{
		double p = urand(rd);
		for (size_t i = 0; i < lights_.size(); ++i)
			if (p < lights_[i].second + eps)
			{
				pmf = lights_[i].second;
				return lights_[i].first;
			}
			else p -= lights_[i].second;

		throw std::logic_error("");
	}

	bool IsVisible (const Vec3 &v, double tol) const
	{ return camera_->IsVisible(v, tol); }

	int Width () const { return width_; }
	int Height () const { return height_; }

	Maybe<Collision> GetCollision (const Ray &r) const
	{
		return collision_->Get(r);
	}
	static Scene Load (const json &json_);
};

#endif //RAYTR_SCENE_H


