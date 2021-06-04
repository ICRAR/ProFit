/**
 * Radial profile base implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Dan Taranu, Rodrigo Tobar
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cmath>
#include <chrono>
#include <map>
#include <sstream>
#include <tuple>
#include <vector>

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/omp_utils.h"
#include "profit/opencl.h"
#include "profit/model.h"
#include "profit/radial.h"
#include "profit/utils.h"

namespace profit
{

inline
void RadialProfile::_image_to_profile_coordinates(double x, double y, double &x_prof, double &y_prof) {
	x -= this->_xcen;
	y -= this->_ycen;
	x_prof =  x * this->_cos_ang + y * this->_sin_ang;
	y_prof = -x * this->_sin_ang_over_axrat + y * this->_cos_ang_over_axrat;
}

double RadialProfile::subsample_pixel(double x0, double x1, double y0, double y1,
                                      unsigned int recur_level, unsigned int max_recursions,
                                      unsigned int resolution) {

	using std::abs;

	double xbin = (x1-x0) / resolution;
	double ybin = (y1-y0) / resolution;
	double half_xbin = xbin/2.;
	double half_ybin = ybin/2.;
	double total = 0;
	double x_prof;
	double y_prof;

	bool recurse = resolution > 1 && recur_level < max_recursions;

#ifdef PROFIT_DEBUG
	/* record how many sub-integrations we've done */
	if( n_integrations.find(recur_level) != n_integrations.end() ) {
		n_integrations[recur_level] += 1;
	}
	else {
		n_integrations[recur_level] = 1;
	}
#endif

	/* The middle X/Y value is used for each pixel */
	std::vector<std::tuple<double, double>> subsample_points;
	double x = x0;

	std::vector<unsigned int> idxs(resolution * resolution);
	if( recurse ) {
		for (unsigned int i = 0; i < resolution; i++) {
			x += half_xbin;
			double y = y0;
			for (unsigned int j = 0; j < resolution; j++) {
				y += half_ybin;

				this->_image_to_profile_coordinates(x, y, x_prof, y_prof);
				double subval = this->evaluate_at(x_prof, y_prof);

				double delta_y_prof = -xbin * _sin_ang_over_axrat + ybin * _cos_ang_over_axrat;
				double testval = this->evaluate_at(abs(x_prof), abs(y_prof) + abs(delta_y_prof));
				if (abs(testval - subval) > acc * subval) {
					subsample_points.emplace_back(std::make_tuple(x, y));
				}
				else {
					total += subval;
				}
				y += half_ybin;
			}

			x += half_xbin;
		}
	}
	else {
		for (unsigned int i = 0; i < resolution; i++) {
			x += half_xbin;
			double y = y0;
			for (unsigned int j = 0; j < resolution; j++) {
				y += half_ybin;
				this->_image_to_profile_coordinates(x, y, x_prof, y_prof);
				total += this->evaluate_at(x_prof, y_prof);
				y += half_ybin;
			}
			x += half_xbin;
		}
	}

	for(auto &point: subsample_points) {
		double x = std::get<0>(point);
		double y = std::get<1>(point);
		total += this->subsample_pixel(x - half_xbin, x + half_xbin,
		                               y - half_ybin, y + half_ybin,
		                               recur_level + 1, max_recursions,
		                               resolution);
	}

	/* Average and return */
	return total / (resolution * resolution);
}

void RadialProfile::initial_calculations() {

	/*
	 * get_rscale() is implemented by subclasses. It provides the translation
	 * from profile-specific parameters into the common rscale concept used in
	 * this common class.
	 */
	this->rscale = this->get_rscale();

	/*
	 * Calculate the total luminosity used by this profile, used
	 * later to calculate the exact contribution of each pixel.
	 */
	double b2 = this->box + 2;
	double r_box = M_PI * b2 / (2 * beta(1/b2, 1/b2));
	double lumtot = this->get_lumtot() * axrat / r_box;
	this->_ie = std::pow(10, -0.4*(this->mag - magzero))/lumtot;

	/*
	 * Optionally adjust the user-given rscale_switch and resolution parameters
	 * to more sensible values that will result in faster profile calculations.
	 */
	if( this->adjust ) {

		/*
		 * Automatially adjust the rscale_switch.
		 * Different profiles do it in different ways
		 */
		this->rscale_switch = this->adjust_rscale_switch();

		/*
		 * Calculate a bound, adaptive upscale
		 */
		auto res = static_cast<unsigned int>(std::ceil(160 / (rscale_switch * rscale)));
		res += res % 2;
		res = std::max(4U, std::min(16U, res));
		this->resolution = res;

		/*
		 * If the user didn't give a rscale_max we calculate one that covers
		 * %99.99 of the flux
		 */
		if( this->rscale_max == 0 ) {
			this->rscale_max = this->adjust_rscale_max();
		}

		/* Adjust the accuracy we'll use for sub-pixel integration */
		acc = adjust_acc(acc);

	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into profile coordinates.
	 *
	 * In galfit the angle started from the Y image axis.
	 */
	double angrad = std::fmod(this->ang + 90, 360.) * M_PI / 180.;
	this->_cos_ang = std::cos(angrad);
	this->_sin_ang = std::sin(angrad);
	_cos_ang_over_axrat = _cos_ang / axrat;
	_sin_ang_over_axrat = _sin_ang / axrat;
}

/**
 * The profile validation function
 */
void RadialProfile::validate() {
	if ( axrat <= 0 ) {
		throw invalid_parameter("axrat <= 0, must have axrat > 0");
	}
	if ( axrat > 1 ) {
		throw invalid_parameter("axrat > 1, must have axrat <= 1");
	}
	if ( box <= -2 ) {
		throw invalid_parameter("box <= -2, must have box > -2");
	}
}

/**
 * The scale by which each image pixel value is multiplied
 */
double RadialProfile::get_pixel_scale(const PixelScale &scale) {
	double pixel_area = scale.x * scale.y;
	return pixel_area * this->_ie;
}

double RadialProfile::adjust_acc(double acc)
{
	return acc;
}

void RadialProfile::subsampling_params(double  /*x*/, double  /*y*/,
                                       unsigned int &resolution,
                                       unsigned int &max_recursions) {
	resolution = this->resolution;
	max_recursions = this->max_recursions;
}

/**
 * The main profile evaluation function
 */
void RadialProfile::evaluate(Image &image, const Mask &mask, const PixelScale &scale,
    const Point &offset, double magzero)
{
	this->magzero = magzero;

	/*
	 * Perform all the pre-calculations needed by the radial profiles
	 * (e.g., Ie, cos/sin ang, etc).
	 * We store these profile-global results in the profile object itself
	 * (it contains extra members to store these values) to avoid passing a long
	 * list of values around every method call.
	 */
	this->initial_calculations();

	// Adjust the center of our profile for the given offset of the image origin
	_xcen = xcen + offset.x * scale.x;
	_ycen = ycen + offset.y * scale.y;

	stats = std::make_shared<RadialProfileStats>();
#ifdef PROFIT_DEBUG
	n_integrations.clear();
#endif /* PROFIT_DEBUG */

#ifndef PROFIT_OPENCL
	evaluate_cpu(image, mask, scale);
#else
	/*
	 * We fallback to the CPU implementation if no OpenCL context has been
	 * given, or if there is no OpenCL kernel implementing the profile
	 */
	auto env = OpenCLEnvImpl::fromOpenCLEnvPtr(model.get_opencl_env());
	if( force_cpu || !env || !supports_opencl() ) {
		evaluate_cpu(image, mask, scale);
		return;
	}

	try {
		if( env->use_double ) {
			evaluate_opencl<double>(image, mask, scale, env);
		}
		else {
			evaluate_opencl<float>(image, mask, scale, env);
		}
	} catch (const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
#endif /* PROFIT_OPENCL */

}

void RadialProfile::evaluate_cpu(Image &image, const Mask &mask, const PixelScale &scale)
{
	double half_xbin = scale.x / 2.;
	double half_ybin = scale.y / 2.;

	auto width = image.getWidth();
	auto height = image.getHeight();
	double flux_scale = this->get_pixel_scale(scale);

	/*
	 * Evaluate the profile at each pixel independently
	 */
	omp_2d_for(model.omp_threads, width, height, [&](unsigned int i, unsigned int j) {

		/* We were instructed to ignore this pixel */
		if( mask && !mask[i + j * width] ) {
			return;
		}

		double x_prof;
		double y_prof;
		double r_prof;
		double y = half_ybin + j * scale.x;
		double x = half_xbin + i * scale.y;
		this->_image_to_profile_coordinates(x, y, x_prof, y_prof);

		/*
		 * Check whether we need further refinement.
		 * TODO: the radius calculation doesn't take into account boxing
		 */
		r_prof = std::sqrt(x_prof*x_prof + y_prof*y_prof);
		double pixel_val;
		if (rscale_max > 0 && r_prof > rscale_max * rscale) {
			pixel_val = 0.;
		}
		else if (rough || r_prof > rscale_switch * rscale) {
			pixel_val = this->evaluate_at(x_prof, y_prof);
		}
		else {

			unsigned int ss_resolution;
			unsigned int ss_max_recursions;
			this->subsampling_params(x, y, ss_resolution, ss_max_recursions);

			/* Subsample and integrate */
			pixel_val =  this->subsample_pixel(x - half_xbin, x + half_xbin,
			                                   y - half_ybin, y + half_ybin,
			                                   0, ss_max_recursions, ss_resolution);
		}

		image[i + j * width] += flux_scale * pixel_val;
	});

}

#ifdef PROFIT_OPENCL

/*
 * Simple structure holding a 2D point
 */
template <typename FT>
struct point_t {
	FT x;
	FT y;
};

/*
 * A structure to hold the information needed to perform subsampling
 * on a specific point
 */
template <typename FT>
class ss_info_t {
public:
	ss_info_t() {};
	ss_info_t(point_t<FT> point, FT xbin, FT ybin, unsigned int resolution, unsigned int max_recursion) :
		point(point), xbin(xbin), ybin(ybin), resolution(resolution), max_recursion(max_recursion) {};
	point_t<FT> point;
	FT xbin;
	FT ybin;
	unsigned int resolution;
	unsigned int max_recursion;
};

/*
 * A similar but smaller version of the ss_info_t structure.
 * This is used as input and output of the subsampling kernel.
 */
template <typename FT>
class ss_kinfo_t {
public:
	ss_kinfo_t() {};
	point_t<FT> point;
	FT xbin;
	FT ybin;
	FT val;
};

template <typename FT>
static inline
unsigned int new_subsampling_points(const std::vector<ss_info_t<FT>> &prev_ss_info, std::vector<ss_info_t<FT>> &ss_info, unsigned int recur_level) {

	ss_info.clear();

	unsigned int subsampled_pixels = 0;
	for(const auto &info: prev_ss_info) {

		const unsigned int res = info.resolution;
		const unsigned int maxr = info.max_recursion;
		FT x = info.point.x;
		FT y = info.point.y;

		if( x == -1 || recur_level > maxr) {
			continue;
		}

		// New subsampling for this point starts at x0,y0
		// and contains res*res subsampling points
		FT x0 = x - info.xbin / 2;
		FT y0 = y - info.ybin / 2;
		FT ss_xbin = info.xbin / FT(res);
		FT ss_ybin = info.ybin / FT(res);

		// we can't cope with more subsampling, sorry
		if( ss_xbin == 0 || ss_ybin == 0 ) {
			continue;
		}

		subsampled_pixels++;
		for(unsigned int j=0; j!=res; j++) {
			FT y_diff = (FT(j) + FT(0.5)) * ss_ybin;
			for(unsigned int i=0; i!=res; i++) {
				FT x_diff = (FT(i) + FT(0.5)) * ss_xbin;
				ss_info.push_back({
					{x0 + x_diff, y0 + y_diff},
					ss_xbin, ss_ybin,
					res, maxr
				});
			}
		}
	}

	return subsampled_pixels;
}

static inline
std::chrono::nanoseconds::rep to_nsecs(const std::chrono::system_clock::duration &d) {
	return std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
}

template <typename FT>
void RadialProfile::evaluate_opencl(Image &image, const Mask & /*mask*/, const PixelScale &scale, OpenCLEnvImplPtr &env) {

	using std::chrono::system_clock;
	typedef point_t<FT> point_t;
	typedef ss_info_t<FT> ss_info_t;
	typedef ss_kinfo_t<FT> ss_kinfo_t;

	unsigned int imsize = image.size();

	OpenCL_times cl_times0 {};
	OpenCL_times ss_cl_times {};
	auto* stats = static_cast<RadialProfileStats *>(this->stats.get());

	/* Points in time we want to measure */
	system_clock::time_point t0, t_kprep, t_opencl, t_loopstart, t_loopend, t_imgtrans;

	/* Prepare the initial evaluation kernel */
	t0 = system_clock::now();
	unsigned int arg = 0;
	auto kname = name + "_" + float_traits<FT>::name;
	cl::Buffer image_buffer = env->get_buffer<FT>(CL_MEM_WRITE_ONLY, imsize);
	cl::Buffer subsampling_points_buffer = env->get_buffer<point_t>(CL_MEM_WRITE_ONLY, imsize);
	cl::Kernel kernel = env->get_kernel(kname);
	kernel.setArg((arg++), image_buffer);
	kernel.setArg((arg++), subsampling_points_buffer);
	kernel.setArg((arg++), image.getWidth());
	kernel.setArg((arg++), image.getHeight());
	kernel.setArg((arg++), (int)rough);
	kernel.setArg((arg++), FT(scale.x));
	kernel.setArg((arg++), FT(scale.y));
	add_common_kernel_parameters<FT>(arg, kernel);
	t_kprep = system_clock::now();

	cl::Event fill_im_evt;
	cl::Event fill_ss_points_evt;
	cl::Event read_evt;

	// OpenCL 1.2 allows to do this; otherwise the work has to be done in the kernel
	// (which we do)
	cl::vector<cl::Event> k_wait_evts;
#if CL_HPP_TARGET_OPENCL_VERSION >= 120
	if( env->get_version() >= 120 ) {
		fill_im_evt = env->queue_fill<FT>(image_buffer, 0);
		fill_ss_points_evt = env->queue_fill<point_t>(subsampling_points_buffer, {-1, -1});
		k_wait_evts.push_back(fill_im_evt);
		k_wait_evts.push_back(fill_ss_points_evt);
	}
#endif /* CL_HPP_TARGET_OPENCL_VERSION >= 120 */

	// Enqueue the kernel, and read back the resulting image + set of points to subsample
	auto kernel_evt = env->queue_kernel(kernel, cl::NDRange(imsize), &k_wait_evts);

	// If FT is double we directly store the result in the profile image
	// Otherwise we have to copy element by element to convert from float to double
	cl::vector<cl::Event> read_waiting_evts{kernel_evt};
	std::vector<FT> image_from_kernel(image.size());
	read_evt = env->queue_read(image_buffer, image_from_kernel.data(), &read_waiting_evts);
	read_evt.wait();
	t_opencl = system_clock::now();
	FT flux_scale = FT(this->get_pixel_scale(scale));
	std::transform(image.begin(), image.end(), image_from_kernel.begin(), image.begin(), [flux_scale](double ipixel, double kpixel) {
		return ipixel + kpixel * flux_scale;
	});
	stats->final_image += to_nsecs(system_clock::now() - t_opencl);

	/* These are the OpenCL-related timings so far */
	cl_times0.kernel_prep = to_nsecs(t_kprep - t0);
	cl_times0.total = to_nsecs(t_opencl - t_kprep);
	if( env->use_profiling ) {
#if CL_HPP_TARGET_OPENCL_VERSION >= 120
		if( env->get_version() >= 120 ) {
			cl_times0.filling_times += cl_cmd_times(fill_im_evt) + cl_cmd_times(fill_ss_points_evt);
		}
#endif /* CL_HPP_TARGET_OPENCL_VERSION >= 120 */
		cl_times0.kernel_times += cl_cmd_times(kernel_evt);
		cl_times0.reading_times += cl_cmd_times(read_evt);
		cl_times0.nwork_items = imsize;
	}

	// we're done here, record the timings and go
	if( rough ) {
		stats->cl_times = std::move(cl_times0);
		stats->total = to_nsecs(system_clock::now() - t0);
		return;
	}

	std::vector<point_t> ss_points(image.size());
	cl::Event read_ss_points_evt = env->queue_read(subsampling_points_buffer, ss_points.data(), &read_waiting_evts);
	read_ss_points_evt.wait();
	if( env->use_profiling ) {
		cl_times0.reading_times += cl_cmd_times(read_ss_points_evt);
	}

	// enrich the points to subsample with their subsampling information
	std::vector<ss_info_t> last_ss_info;
	last_ss_info.reserve(image.size());
	unsigned int top_recursions = 0;
	for(auto const &point: ss_points) {
		if( point.x == -1 ) {
			continue;
		}
		unsigned int resolution;
		unsigned int max_recursions;
		subsampling_params(point.x, point.y, resolution, max_recursions);
		top_recursions = std::max(top_recursions, max_recursions);
		last_ss_info.push_back({point, FT(scale.x), FT(scale.y), resolution, max_recursions});
	}

	auto ss_kname = name + "_subsample_" + float_traits<FT>::name;
	cl::Kernel subsample_kernel = env->get_kernel(ss_kname);

	class im_result_t {
	public:
		im_result_t() {};
		point_t point;
		FT value;
	};
	std::vector<im_result_t> subimages_results;

	// Preparing for the recursive subsampling
	std::vector<ss_kinfo_t> ss_kinfo;
	std::vector<ss_info_t> ss_info;
	unsigned int recur_level = 0;

	t_loopstart = system_clock::now();
	unsigned int subimage_idx = 0;
	while( recur_level <= top_recursions ) {

		/* Points in time we want to measure */
		system_clock::time_point t0, t_newsamples, t_trans_h2k, t_kprep, t_opencl, t_trans_k2h;


		t0 = system_clock::now();
		unsigned int subsampled_pixels = new_subsampling_points<FT>(last_ss_info, ss_info, recur_level);
		t_newsamples = system_clock::now();

		auto subsamples = ss_info.size();
		if( !subsamples ) {
			break;
		}

		ss_cl_times.nwork_items += subsamples;
#ifdef PROFIT_DEBUG
		/* record how many sub-integrations we've done */
		n_integrations[recur_level] = subsampled_pixels;
#else
		// avoid warning because of unused variable
		UNUSED(subsampled_pixels);
#endif

		/* Keeping things in size */
		auto prev_im_size = subimages_results.size();
		subimages_results.resize(prev_im_size + subsamples);
		last_ss_info.resize(subsamples);

		try {

			cl::Buffer ss_kinfo_buf = env->get_buffer<ss_kinfo_t>(CL_MEM_READ_WRITE, subsamples);

			arg = 0;
			subsample_kernel.setArg((arg++), ss_kinfo_buf);
			subsample_kernel.setArg((arg++), FT(acc));
			add_common_kernel_parameters<FT>(arg, subsample_kernel);

			t_kprep = system_clock::now();

			// The information we pass down to the kernels is a subset of the original
			ss_kinfo.resize(subsamples);
			std::transform(ss_info.begin(), ss_info.end(), ss_kinfo.begin(), ss_kinfo.begin(), [](const ss_info_t &info, ss_kinfo_t &kinfo) -> ss_kinfo_t {
				kinfo.point = info.point;
				kinfo.xbin = info.xbin;
				kinfo.ybin = info.ybin;
				return kinfo;
			});
			t_trans_h2k = system_clock::now();

			auto w_ss_kinfo_evt = env->queue_write(ss_kinfo_buf, ss_kinfo.data());
			cl::vector<cl::Event> kernel_waiting_evts{w_ss_kinfo_evt};
			auto ss_kernel_evt = env->queue_kernel(subsample_kernel, cl::NDRange(subsamples), &kernel_waiting_evts);
			cl::vector<cl::Event> read_waiting_evts{ss_kernel_evt};
			auto r_ss_kinfo_evt = env->queue_read(ss_kinfo_buf, ss_kinfo.data(), &read_waiting_evts);
			r_ss_kinfo_evt.wait();
			t_opencl = system_clock::now();

			// Feed back the kinfo to the main subsampling info vectors
			auto ss_info_it = ss_info.begin();
			auto last_ss_info_it  = last_ss_info.begin();
			for(auto &kinfo: ss_kinfo) {

				// Copy the point information from the kernel
				last_ss_info_it->point = kinfo.point;

				// ... and the rest of the subsampling info for next round
				last_ss_info_it->xbin = ss_info_it->xbin;
				last_ss_info_it->ybin = ss_info_it->ybin;
				last_ss_info_it->resolution = ss_info_it->resolution;
				last_ss_info_it->max_recursion = ss_info_it->max_recursion;

				FT val = kinfo.val;
				for(unsigned int i=0; i<=recur_level; i++) {
					val /= FT(ss_info_it->resolution * ss_info_it->resolution);
				}

				subimages_results[subimage_idx].point = std::move(ss_info_it->point);
				subimages_results[subimage_idx].value = val;
				subimage_idx++;

				last_ss_info_it++;
				ss_info_it++;
			}
			t_trans_k2h = system_clock::now();

			stats->subsampling.new_subsampling += to_nsecs(t_newsamples - t0);
			stats->subsampling.inital_transform += to_nsecs(t_trans_h2k - t_kprep);
			stats->subsampling.final_transform += to_nsecs(t_trans_k2h - t_opencl);
			ss_cl_times.kernel_prep += to_nsecs(t_kprep - t_newsamples);
			ss_cl_times.total += to_nsecs(t_opencl - t_trans_h2k);
			if( env->use_profiling ) {
				ss_cl_times.kernel_times += cl_cmd_times(ss_kernel_evt);
				ss_cl_times.writing_times += cl_cmd_times(w_ss_kinfo_evt);
				ss_cl_times.reading_times += cl_cmd_times(r_ss_kinfo_evt);
			}

		} catch(const cl::Error &e) {
			// running out of memory, cannot go any further
			if( e.err() == CL_INVALID_BUFFER_SIZE ) {
				break;
			}
			throw;
		}

		recur_level++;
	}

	t_loopend = system_clock::now();

	std::for_each(subimages_results.begin(), subimages_results.end(), [&image, &scale, &flux_scale](const im_result_t &res) {
		FT x = FT(res.point.x / scale.x);
		FT y = FT(res.point.y / scale.y);
		unsigned int idx = static_cast<unsigned int>(floor(x)) + static_cast<unsigned int>(floor(y)) * image.getWidth();
		image[idx] += res.value * flux_scale;
	});

	// the image needs to be multiplied by the pixel scale
	t_imgtrans = system_clock::now();

	stats->subsampling.pre_subsampling = to_nsecs(t_loopstart - t_opencl);
	stats->subsampling.total = to_nsecs(t_loopend - t_loopstart);
	stats->final_image += to_nsecs(t_imgtrans - t_loopend);
	stats->total = to_nsecs(t_imgtrans - t0);
	stats->cl_times = std::move(cl_times0);
	stats->subsampling.cl_times = std::move(ss_cl_times);

}

template <typename FT>
void RadialProfile::add_common_kernel_parameters(unsigned int arg, cl::Kernel &kernel) const {
	kernel.setArg((arg++), FT(_xcen));
	kernel.setArg((arg++), FT(_ycen));
	kernel.setArg((arg++), FT(_cos_ang));
	kernel.setArg((arg++), FT(_sin_ang));
	kernel.setArg((arg++), FT(axrat));
	kernel.setArg((arg++), FT(rscale));
	kernel.setArg((arg++), FT(rscale_switch));
	kernel.setArg((arg++), FT(rscale_max));
	kernel.setArg((arg++), FT(box));
	if( float_traits<FT>::is_float ) {
		add_kernel_parameters_float(arg, kernel);
	}
	else {
		add_kernel_parameters_double(arg, kernel);
	}
}

bool RadialProfile::supports_opencl() const {
	return true;
}

#endif /* PROFIT_OPENCL */

/**
 * Constructor with sane defaults
 */
RadialProfile::RadialProfile(const Model &model, const std::string &name) :
	Profile(model, name),
	xcen(0), ycen(0),
	mag(15), ang(0),
	axrat(1), box(0),
	rough(false), acc(0.1),
	rscale_switch(1), resolution(9),
	max_recursions(2), adjust(true),
	rscale_max(0),
	force_cpu(false),
	rscale(0), _ie(0),
	_cos_ang(0), _sin_ang(0), _cos_ang_over_axrat(0), _sin_ang_over_axrat(0),
	magzero(0)
{
	register_parameter("rough", rough);
	register_parameter("adjust", adjust);
	register_parameter("force_cpu", force_cpu);
	register_parameter("xcen", xcen);
	register_parameter("ycen", ycen);
	register_parameter("mag", mag);
	register_parameter("ang", ang);
	register_parameter("axrat", axrat);
	register_parameter("box", box);
	register_parameter("acc", acc);
	register_parameter("rscale_switch", rscale_switch);
	register_parameter("rscale_max", rscale_max);
	register_parameter("max_recursions", max_recursions);
	register_parameter("resolution", resolution);
}

#ifdef PROFIT_DEBUG
std::map<int,int> RadialProfile::get_integrations() {
	return n_integrations;
}
#endif

#ifdef PROFIT_OPENCL
void RadialProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	// subclasses may add more kernel parameters
}

void RadialProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	// subclasses may add more kernel parameters
}
#endif /* PROFIT_OPENCL */

} /* namespace profit */