//
// Created by dc on 5/17/16.
//

#include "objekt.h"
#include "shade.h"

Primitive::Primitive (PSurfaceProperty surfaceProperty):
		surfaceProperty_(surfaceProperty), light_property_(nullptr)
{
	is_light_ = surfaceProperty_->emission.nonzero();
	if (is_light_)
		SetLightProperty(make_shared<LightProperty>());
}

void Primitive::SetLightProperty (PLightProperty prop)
{
	light_property_ = prop;
	light_property_->prim_ = this;
}

void Primitive::SetUnidirectedLight (Vec3 i, Vec3 j)
{
	light_property_->light_i_ = i;
	light_property_->light_j_ = j;
	light_property_->unidirected_light_ = true;
}

void Primitive::SetLightConcentration (int e)
{
	light_property_->light_exp_ = e;
}

void Primitive::SampleLightDir (rand_data *rd, Vec3 &pos, Vec3 &N, Vec3 &dir, double &pos_invpdf, double &dir_pdf) const
{
	// TODO: refactor
	light_property_->SampleDir(rd, pos, N, dir, pos_invpdf, dir_pdf);
}

