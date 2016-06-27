//
// Created by dc on 4/7/16.
//

#ifndef RAYTR_OBJEKT_H
#define RAYTR_OBJEKT_H

#include "common.h"
#include "random.h"

class Primitive;

struct Collision
{
	double t; // pos = S + tD
	Vec3 pos, norm;
	const Primitive *obj;
	bool inside; // true := diffusion = n_t/n_i

	Collision () {}
	Collision (double t_, Vec3 P, Vec3 N, const Primitive *obj_, bool inside_):
			t(t_), pos(P), norm(N), obj(obj_), inside(inside_) {}
};

struct BBox
{
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	BBox (double xmn = 0., double ymn = 0., double zmn = 0.,
		  double xmx = 0., double ymx = 0., double zmx = 0.):
	xmin(xmn), ymin(ymn), zmin(zmn), xmax(xmx), ymax(ymx), zmax(zmx) {}
	BBox (const Vec3 &vmin, const Vec3 &vmax) {
		xmin = vmin.x;
		ymin = vmin.y;
		zmin = vmin.z;
		xmax = vmax.x;
		ymax = vmax.y;
		zmax = vmax.z;
	}
	static BBox Min () { return BBox(inf, inf, inf, -inf, -inf, -inf); }
	static BBox Max () { return BBox(-inf, -inf, -inf, inf, inf, inf); }
	inline double & crd_min (int p) { return (&xmin)[p]; }
	inline double & crd_max (int p) { return (&xmax)[p]; }
	inline double crd_min (int p) const { return (&xmin)[p]; }
	inline double crd_max (int p) const { return (&xmax)[p]; }
	inline int max_axis () const
	{
		double d[3] = {max(0., xmax - xmin), max(0., ymax - ymin), max(0., zmax - zmin)};
		return int(std::max_element(d, d + 3) - d);
	}
	inline void merge_with (const BBox &b)
	{
		if (b.xmin < xmin) xmin = b.xmin;
		if (b.ymin < ymin) ymin = b.ymin;
		if (b.zmin < zmin) zmin = b.zmin;
		if (b.xmax > xmax) xmax = b.xmax;
		if (b.ymax > ymax) ymax = b.ymax;
		if (b.zmax > zmax) zmax = b.zmax;
	}
	inline BBox merged_with (const BBox &b)
	{
		return BBox{min(xmin, b.xmin), min(ymin, b.ymin), min(zmin, b.zmin),
				    max(xmax, b.xmax), max(ymax, b.ymax), max(zmax, b.zmax)};
	}
	inline double SA ()
	{
		double dx = max(0., xmax - xmin),
			   dy = max(0., ymax - ymin),
			   dz = max(0., zmax - zmin);
		return sqr(dx + dy + dz) - sqr(dx) - sqr(dy) - sqr(dz);
	}
};

struct ShadingProperty;
typedef shared_ptr<const ShadingProperty> PSurfaceProperty;
struct LightProperty;
typedef shared_ptr<LightProperty> PLightProperty;
typedef shared_ptr<const LightProperty> CPLightProperty;

class Primitive
{
protected:
	PSurfaceProperty surfaceProperty_;
	bool is_light_;
	double area_;
	BBox bbox_;
	PLightProperty light_property_;
public:
	Primitive (PSurfaceProperty surfaceProperty);

	// @param dir: unit vector.
	virtual Maybe<Collision> GetCollision (const Vec3 &src, const Vec3 &dir) const = 0;

	// Return the solid angle extended by *this seen from src, which is the
	// inversion of light pdf in direct light calculation.
	virtual double SAInvPDF (const Vec3 &src) const = 0;

	// Sample a point on surface that is not sheltered when viewed from src.
	virtual Vec3 SampleSurface (rand_data *rd, const Vec3 &src) const = 0;
	// Sample pair(surface-point, local-normal).
	virtual void SampleSurface (rand_data *rd, Vec3 &sample, Vec3 &N) const = 0;

	// Convert world coordinate to UV (texture) coordinate.
	// Throws domain_error if unimplemented; check surfaceProperty_ before calling this.
	virtual Vec3 ToUV (const Vec3 &surfacePoint) const = 0;

	virtual double Area () const
	{
		return area_;
	}
	virtual const PSurfaceProperty &GetSurfaceProperty() const
	{
		return surfaceProperty_;
	}
	bool IsLight() const
	{
		return is_light_;
	}
	BBox GetBBox () const
	{
		return bbox_;
	}

	void SetLightProperty (PLightProperty prop);

	CPLightProperty GetLightProperty () const
	{
		return light_property_;
	}
	void SetUnidirectedLight (Vec3 i, Vec3 j);

	void SetLightConcentration (int e);

	void SampleLightDir (rand_data *rd, Vec3 &pos, Vec3 &N, Vec3 &dir, double &pos_invpdf, double &dir_pdf) const;

	// AABB GetAABB() = 0;
};

#endif //RAYTR_OBJEKT_H
