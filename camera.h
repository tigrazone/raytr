//
// Created by dc on 4/7/16.
//

#ifndef RAYTR_CAMERA_H
#define RAYTR_CAMERA_H

#include "common.h"

class Camera
{
public:
	// Cast ray with respect to screen point; return pair<src, normalized_dir>
	virtual void Project (rand_data *rd, double Xscr, double Yscr, Ray &ray) const = 0;

	// Get camera origin for debugging
	virtual Vec3 O () const = 0;

	// Check if a sphere (P, rad) has visible part on screen
	virtual bool IsVisible (Vec3 P, double rad) const = 0;

	// Return the position of P on screen, if visible
	virtual Maybe<Vec3> Rasterize (Vec3 P) const = 0;
};

class PCamera: public Camera
{
	Vec3 O_, Oimg_;                  // camera and its projection on image plane in world coordinate
	Vec3 I_, J_, K_;                 // {x, y, z axis of camera space} in world coordinate
	double camWindow_[4];               // range of x and y on image plane
	double imgW_, imgH_;                // size of image
	double aperture_, focal_distance_;  // for depth of field effect
	bool aa_;

public:
	PCamera (Vec3 O, Vec3 Oimg,
			 Vec3 I, Vec3 J,
			 const double camWindow[4],
			 int imgW, int imgH,
			 double aperture, double fod,
			 bool aa = false):
			O_(O), Oimg_(Oimg),
			I_(I.normalized()), J_(J.normalized()),
			imgW_(imgW), imgH_(imgH),
			aperture_(aperture), focal_distance_(fod), aa_(aa)
	{
		K_ = I_.cross(J_);
		std::copy(camWindow, camWindow + 4, camWindow_);
	}

	PCamera (const PCamera &r)
	{
		memcpy(this, &r, sizeof(PCamera));
	}

	void Project (rand_data *rd, double Xscr, double Yscr, Ray &ray) const override
	{
		Vec3 Pndc(Xscr / imgW_ * (camWindow_[1] - camWindow_[0]) + camWindow_[0],
				  Yscr / imgH_ * (camWindow_[3] - camWindow_[2]) + camWindow_[2],
				  0.);
		Vec3 Pw = Pndc.x * I_ + Pndc.y * J_ + Oimg_;
		ray.S = O_;
		ray.D = (Pw - ray.S).normalized();

		// NOTE: For compatibility, we assume the image plane is between camera 
		// scene if aperture_ < eps (while the opposite for lens effect). So 
		// O and Oimg should be set differently in two cases. 
		// TODO: deprecate it.
		if (aperture_ > eps)
		{
			Vec3 Pfp = O_ + (focal_distance_ / ray.D.dot(K_)) * ray.D;
			double dx, dy;
			URandPlate(rd, aperture_, dx, dy);
			ray.S += dx * I_ + dy * J_;
			ray.D = (Pfp - ray.S).normalized();
		}
	}

	virtual Maybe<Vec3> Rasterize (Vec3 P) const override
	{
		if ((P - O_).dot(Oimg_ - O_) < eps)
			return Maybe<Vec3>::None();

		auto dir = (P - O_).normalized();
		auto z = dir.dot((Oimg_ - O_).normalized());
		auto Pcp = O_ + ((Oimg_ - O_).norm()) / z * dir - Oimg_;
		double x = Pcp.dot(I_), y = Pcp.dot(J_);
		x = int((x - camWindow_[0]) / (camWindow_[1] - camWindow_[0]) * imgW_);
		y = int((y - camWindow_[2]) / (camWindow_[3] - camWindow_[2]) * imgH_);
		if (x < eps || y < eps || x + 1 > imgW_ || y + 1 > imgH_)
			return Maybe<Vec3>::None();

		return Maybe<Vec3>::Some(Vec3(x, y, 0));
	}

	virtual bool IsVisible (Vec3 P, double rad) const override
	{
		// The most probable position is near camera
		return Rasterize(P - rad * (P - O_).normalized()).valid;
	}


	Vec3 O () const override
	{
		return O_;
	}


	static PCamera* Parse (json j, int imgWidth, int imgHeight)
	{
		auto jCamWindow = j.at("cameraWindow");
		double camWindow[4];
		for (int t = 0; t < 4; ++t)
			camWindow[t] = jCamWindow.at(t).get<double>();

		return new PCamera(
				Vec3::Parse(j.at("O")),
				Vec3::Parse(j.at("Oimg")),
				Vec3::Parse(j.at("I")),
				Vec3::Parse(j.at("J")),
				camWindow,
				imgWidth, imgHeight,
				J_AT_F(j, "aperture"),
				J_AT_F(j, "fod"));
	}
};

// TODO: deprecated. remove me.
class LCamera: public Camera
{
	Vec3 o_;
public:
	LCamera (Vec3 o) : o_(o)
	{ }

	LCamera ()
	{ }

	Vec3 O () const override
	{ return o_; }

	void Project (rand_data *rd, double x, double y, Ray &ray) const override
	{
		ray.S = o_;
		ray.D = (Vec3(x, y, 0) - o_).normalized();
	}

	virtual Maybe<Vec3> Rasterize (Vec3 P) const
	{
		throw std::logic_error("TODO");
	}

	virtual bool IsVisible (Vec3 P, double rad) const override
	{
		throw std::logic_error("TODO");
	}

};

#endif //RAYTR_CAMERA_H
