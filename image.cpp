//
// Created by dc on 4/11/16.
//

#include "image.h"

#ifdef NO_CV
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"
#endif

using namespace std;

Image::Image (int w, int h, const char *dump_prefix, int display)
		: w_(w), h_(h), finished(0), dump_prefix_(dump_prefix)
{
	orig_data_.resize(h);
	for (int i = 0; i < h; ++i)
		orig_data_[i] = vector<Spectrum>(w, Spectrum(0.));

#ifndef NO_CV
	data_ = cv::Mat::zeros(h, w, CV_8UC3);
	if (display)
	{
		cv::namedWindow(dump_prefix_.c_str());
		displayer_ = new std::thread([this](){
			while (1)
			{
				cv::imshow(dump_prefix_.c_str(), data_);
				cv::waitKey(200);
				if (this->finished)
					break;
			}
		});
	}
#else
	data_ = new u_char[h * w * 3];
	fill(data_, data_ + h * w * 3, u_char(0));
#endif
}


void Image::Set (int x, int y, const Spectrum &color_)
{
	assert(!(-1 * color_).nonzero());

	RGBSpectrum color = color_.ToRGB();
#ifdef SPECTRAL
	// In spectral mode, if sample size is too small, we may get negative RGB value
	color.r = max(color.r, 0.);
	color.g = max(color.g, 0.);
	color.b = max(color.b, 0.);
#endif

	orig_data_[h_ - 1 - y][x] = color_;

	// Draw on canvas
#define E(x) pow(400.0 * (1.0 - exp(-x / 400.0)), 1.0 / 2.2)
	color.b = E(color.b);
	color.r = E(color.r);
	color.g = E(color.g);

	color.b = std::min(color.b, 1.);
	color.r = std::min(color.r, 1.);
	color.g = std::min(color.g, 1.);
	color = color.operator*(255);
#ifndef NO_CV
	data_.at<cv::Vec3b>(h_ - 1 - y, x) = cv::Vec3b((uchar)color.b, (uchar)color.g, (uchar)color.r);
#else
	int offs = (h_ - 1 - y) * w_ * 3 + x * 3;
	data_[offs + 0] = (u_char)color.r;
	data_[offs + 1] = (u_char)color.g;
	data_[offs + 2] = (u_char)color.b;
#endif
}

void Image::Show ()
{
	finished = 1;
}

static void _Histogram (const vector<vector<Spectrum>> &);

void Image::Dump (int id)
{
	std::ostringstream fileName_;
	fileName_ << dump_prefix_ << "." << id << ".png";
	auto fileName = fileName_.str();
	const char *fileNameC = fileName.c_str();
#ifndef NO_CV
	for (int x = 0; x < w_; x += 60)
	{
		cv::line(data_, cv::Point(x, 0), cv::Point(x, 20), cv::Scalar(200, 200, 200, 240));
	}
	for (int y = 0; y < h_; y += 60)
	{
		cv::line(data_, cv::Point(0, y), cv::Point(20, y), cv::Scalar(200, 200, 200, 240));
	}
	cv::imwrite(fileNameC, data_);
#else
	stbi_write_png(fileNameC, w_, h_, 3, (void*)data_, 0);
#endif
	clog << "Image dumped to " << fileName << ".\n";
	// Do some stats.
	_Histogram(orig_data_);
}

static void _Histogram (const vector<vector<Spectrum>> &orig_data_)
{
#ifdef NO_CV
	return;
#endif

	int h_ = (int)orig_data_.size();
	int w_ = (int)orig_data_[0].size();
	int n_bins = 10, n_elems = w_ * h_ / 2;
	int hhgt = 400, hwid = 40;

	double max_co = 0.;
	for (int d = 0; d < Spectrum::N; ++d)
	{
		for (int x = 0; x < w_ / 2; ++x)
			for (int y = 0; y < h_; ++y)
			{
				double cr = orig_data_[y][x][d];
				if (cr > max_co)
					max_co = cr;
			}
	}
	clog << "\tMax intensity = " << max_co << "\n";

	double histogram[Spectrum::N][n_bins], hist_ruler[n_bins];
	memset(histogram, 0, sizeof(histogram));

	hist_ruler[n_bins - 1] = 100.0;
	for (int t = n_bins - 2; t >= 0; --t)
		hist_ruler[t] = hist_ruler[t + 1] * 0.25;

	clog << "\tHistogram ruler = 0 " << std::scientific << std::setprecision(2);
	for (int t = 0; t < n_bins; ++t)
		clog << hist_ruler[t] << " \n"[t + 1 == n_bins];

	for (int ch = 0; ch < Spectrum::N; ++ch)
	{
		for (int x = 0; x < w_ / 2; ++x)
			for (int y = 0; y < h_; ++y)
			{
				double cr = orig_data_[y][x][ch];
				int k;
				for (k = 0; k + 1 < n_bins; ++k)
					if (cr < hist_ruler[k])
						break;
				histogram[ch][k] += 1.0 / n_elems;
			}
	}

	clog << std::fixed;
	for (int d = 0; d < Spectrum::N; ++d)
	{
		clog << "\tHist " << d << ": ";
		for (int k = 0; k < n_bins; ++k)
			clog << histogram[d][k] * 1e5 << " \n"[k + 1 == n_bins];
	}

#ifdef DISP_HISTOGRAM
	static bool first_run_ = true;

	if (first_run_)
	{
		cv::namedWindow("Histogram");
		first_run_ = false;
	}

	cv::Mat image;
	image = cv::Mat::zeros(hhgt, hwid * n_bins, CV_8UC3);
	for (int d = 0; d < 3; ++d)
	{
		for (int k = 0; k + 1 < n_bins; ++k)
		{
			cv::line(image,
					 cv::Point(k * hwid, histogram[d][k] * hhgt),
					 cv::Point((k + 1) * hwid, histogram[d][k + 1] * hhgt),
					 cv::Scalar(200 * (d == 0), 200 * (d == 1), 200 * (d == 2), 200));
		}
	}
	cv::imshow("Histogram", image);
	cv::waitKey(-1);
#endif
}

