//
// Created by dc on 4/7/16.
//

#include "scene.h"
#include "primitives.h"
#include "bbox_tree.h"
#include "objloader.h"

Scene Scene::Load (const json &json_)
{
	if (json_.count("lights"))
		throw std::domain_error("lights field is deprecated");

	int width = json_.at("w").get<int>(),
		height = json_.at("h").get<int>();

	// Create camera
	auto jCamera = json_.at("camera");

	Camera *camera_;
	if (!jCamera.is_array())
		camera_ = PCamera::Parse(jCamera, width, height);
	else 
		camera_ = new LCamera(Vec3::Parse(json_.at("camera")));

	// Create objects
	std::vector<Primitive*> lights_;
	std::vector<Primitive*> objects_;
	std::vector<Destroyable*> resources_;

	if (json_.count("obj_model"))
	{
		int interpolate_normal = json_.count("interpolate_normal") ?
								  json_.at("interpolate_normal").get<int>() :
								  1;
		json jmodifiers = json_.count("obj_mod") ? json_["obj_mod"] : json();
		LoadMesh(json_.at("obj_model").get<string>(),
				 jmodifiers,
				 interpolate_normal,
				 objects_,
				 resources_);
	}
	if (json_.count("objects"))
	{
		auto j_objects = json_.at("objects");
		for (const json &jo: j_objects)
		{
			Primitive *object = CreateBasicObject(jo);
			objects_.push_back(object);
		}
	}
	for (auto object: objects_)
	{
		if (object->IsLight())
			lights_.push_back(object);
	}

	Spectrum env_light(0);
	if (json_.count("env_light"))
		env_light = Spectrum::Parse(json_["env_light"]);

	auto colMgr = new BBoxTree(objects_);
//	auto colMgr = new BBTChecker(objects_);
//	auto colMgr = new BFCollisionMGR(objects_);
	return Scene(objects_, lights_,
				 resources_,
				 colMgr, camera_,
				 env_light, width, height);
}

PointLight *PointLight::Create (const json &j)
{
	auto O = Vec3::Parse(j.at("O"));
	auto intensity = Spectrum::Parse(j.at("I"));
	return new PointLight(O, intensity);
}


