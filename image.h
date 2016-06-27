//
// Created by dc on 4/8/16.
//

#ifndef RAYTR_IMAGE_H
#define RAYTR_IMAGE_H

#include "common.h"
#include "spectrum.h"

#include <thread>

// #define NO_CV

#ifndef NO_CV
#include <opencv2/opencv.hpp>
#else
// Use stb_image_write; included in .cpp.
#endif

struct Image
{
private:
	int w_, h_, finished;
#ifndef NO_CV
	cv::Mat data_;
#else
	u_char* data_;
#endif
	std::vector<std::vector<Spectrum>> orig_data_; // May be >1; use to draw histogram
	std::thread *displayer_;
	std::string dump_prefix_;
public:
	Image (int w, int h, const char *dump_prefix, int display);
	void Set (int x, int y, const Spectrum &color);
	void Show ();
	void Dump (int id);
};

#endif //RAYTR_IMAGE_H
