/**
 * Image convolution implementation
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
#include <functional>
#include <memory>
#include <sstream>
#include <vector>

#include "profit/convolve.h"
#include "profit/exceptions.h"
#include "profit/utils.h"


namespace profit
{

Convolver::~Convolver()
{
	// no-op
}


Image BruteForceConvolver::convolve(const Image &src, const Image &krn, const Mask &mask)
{

	auto src_width = src.getWidth();
	auto src_height = src.getHeight();
	auto krn_width = krn.getWidth();
	auto krn_height = krn.getHeight();

	double pixel;
	unsigned int i, j, k, l;
	unsigned int krn_half_width = krn_width / 2;
	unsigned int krn_half_height = krn_height / 2;
	unsigned int krn_size = krn_width * krn_height;
	int src_i, src_j;

	Image convolution(src_width, src_height);

	const double *krn_data = krn.getData().data();
	double *out = convolution.getData().data() - 1;
	const double *srcPtr1 = src.getData().data() - 1, *srcPtr2;
	const double *krnPtr;
	auto mask_it = mask.getData().begin();

	/* Convolve! */
	/* Loop around the output image first... */
	for (j = 0; j < src_height; j++) {
		for (i = 0; i < src_width; i++) {

			out++;
			srcPtr1++;

			/* Don't convolve this pixel */
			if( !mask.empty() ) {
				if( !*mask_it++ ) {
					*out = 0;
					continue;
				}
			}

			pixel = 0;
			krnPtr = krn_data + krn_size - 1;
			srcPtr2 = srcPtr1 - krn_half_width - krn_half_height*src_width;

			/* ... now loop around the kernel */
			for (l = 0; l < krn_height; l++) {

				src_j = (int)j + (int)l - (int)krn_half_height;
				for (k = 0; k < krn_width; k++) {

					src_i = (int)i + (int)k - (int)krn_half_width;

					if( src_i >= 0 && (unsigned int)src_i < src_width &&
					    src_j >= 0 && (unsigned int)src_j < src_height ) {
						pixel +=  *srcPtr2 * *krnPtr;
					}

					srcPtr2++;
					krnPtr--;
				}
				srcPtr2 += src_width - krn_width;
			}

			*out = pixel;
		}
	}

	return convolution;
}

#ifdef PROFIT_FFTW
FFTConvolver::FFTConvolver(unsigned int src_width, unsigned int src_height,
                           unsigned int krn_width, unsigned int krn_height,
                           FFTPlan::effort_t effort, unsigned int plan_omp_threads,
                           bool reuse_krn_fft) :
	plan(),
	krn_fft(),
	reuse_krn_fft(reuse_krn_fft)
{

	if (krn_width > src_width) {
		throw invalid_parameter("krn_width must be <= src_width");
	}
	if (krn_height > src_height) {
		throw invalid_parameter("krn_height must be <= src_height");
	}
	auto convolution_size = 4 * src_width * src_height;
	plan = std::unique_ptr<FFTPlan>(new FFTPlan(convolution_size, effort, plan_omp_threads));
}

Image FFTConvolver::convolve(const Image &src, const Image &krn, const Mask &mask)
{

	typedef std::complex<double> complex;

	auto src_width = src.getWidth();
	auto src_height = src.getHeight();
	auto krn_width = krn.getWidth();
	auto krn_height = krn.getHeight();

	// Create extended images first
	auto ext_width = 2 * src_width;
	auto ext_height = 2 * src_height;
	Image ext_img = src.extend(ext_width, ext_height, 0, 0);

	// Forward FFTs
	std::vector<complex> src_fft = plan->forward(ext_img);
	if (krn_fft.empty()) {
		auto krn_start_x = (src_width - krn_width) / 2;
		auto krn_start_y = (src_height - krn_height) / 2;
		Image ext_krn = krn.extend(ext_width, ext_height, krn_start_x, krn_start_y);
		krn_fft = plan->forward(ext_krn);
	}

	// element-wise multiplication
	std::transform(src_fft.begin(), src_fft.end(), krn_fft.begin(), src_fft.begin(), std::multiplies<complex>());

	if (!reuse_krn_fft) {
		krn_fft.clear();
	}

	// inverse FFT and scale down
	Image res(plan->backward_real(src_fft), ext_width, ext_height);
	res /= res.getSize();

	// crop the image to original size, apply the mask, and good bye
	auto x_offset = src_width / 2;
	auto y_offset = src_height / 2;

	// even image and odd kernel requires slight adjustment
	if (src_width % 2 == 0 or krn_width % 2 == 0) {
		x_offset -= 1;
	}
	if (src_height % 2 == 0 or krn_height % 2 == 0) {
		y_offset -= 1;
	}

	auto cropped = res.crop(src_width, src_height, x_offset, y_offset);
	cropped &= mask;
	return cropped;
}

#endif /* PROFIT_FFTW */

#ifdef PROFIT_OPENCL
OpenCLConvolver::OpenCLConvolver(std::shared_ptr<OpenCL_env> opencl_env) :
	env(opencl_env)
{
	if (!env) {
		throw invalid_parameter("Empty OpenCL environment given to OpenCLConvolver");
	}
}

Image OpenCLConvolver::convolve(const Image &src, const Image &krn, const Mask &mask)
{
	try {
		return _convolve(src, krn, mask);
	} catch (const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error while convolving: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

Image OpenCLConvolver::_convolve(const Image &src, const Image &krn, const Mask &mask) {

	// We use a group size of 16x16, so let's extend the src image
	// to the next multiple of 16
	auto clpad_x = (16 - (src.getWidth() % 16)) % 16;
	auto clpad_y = (16 - (src.getHeight() % 16)) % 16;
	const Image clpad_src = src.extend(src.getWidth() + clpad_x,
	                                   src.getHeight() + clpad_y, 0, 0);

	// Convolve using the appropriate data type
	Image result;
	if (env->use_double) {
		result = _clpadded_convolve<double>(clpad_src, krn, src);
	}
	else {
		result = _clpadded_convolve<float>(clpad_src, krn, src);
	}

	// Crop the resulting image, mask
	return result.crop(src.getWidth(), src.getHeight(), 0, 0) & mask;
}

template<typename T>
Image OpenCLConvolver::_clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src) {

	using cl::Buffer;
	using cl::Event;
	using cl::Kernel;
	using cl::NDRange;
	using cl::NullRange;

	Buffer src_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, src.getSize());
	Buffer krn_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, krn.getSize());
	Buffer conv_buf = env->get_buffer<T>(CL_MEM_WRITE_ONLY, src.getSize());

	std::vector<T> src_data(src.getSize());
	std::copy(src.getData().begin(), src.getData().end(), src_data.begin());
	std::vector<T> krn_data(krn.getSize());
	std::copy(krn.getData().begin(), krn.getData().end(), krn_data.begin());

	// Write both images' data to the device
	Event src_wevt = env->queue_write(src_buf, src_data.data());
	Event krn_wevt = env->queue_write(krn_buf, krn_data.data());

	// We need this much local memory on each local group
	auto local_size = sizeof(T);
	local_size *= (16 + 2 * (krn.getWidth() / 2));
	local_size *= (16 + 2 * (krn.getHeight() / 2));

	// Prepare the kernel
	auto kname = std::string("convolve_") + float_traits<T>::name;
	Kernel clKernel = env->get_kernel(kname);
	clKernel.setArg(0, src_buf);
	clKernel.setArg(1, orig_src.getWidth());
	clKernel.setArg(2, orig_src.getHeight());
	clKernel.setArg(3, krn_buf);
	clKernel.setArg(4, krn.getWidth());
	clKernel.setArg(5, krn.getHeight());
	clKernel.setArg(6, conv_buf);

	// Execute
	std::vector<Event> exec_wait_evts {src_wevt, krn_wevt};
	auto exec_evt = env->queue_kernel(clKernel, NDRange(src.getWidth(), src.getHeight()), &exec_wait_evts, NDRange(16, 16));

	// Read and good bye
	std::vector<Event> read_wait_evts {exec_evt};
	std::vector<T> conv_data(src.getSize());
	Event read_evt = env->queue_read(conv_buf, conv_data.data(), &read_wait_evts);
	read_evt.wait();

	Image conv(src.getWidth(), src.getHeight());
	std::copy(conv_data.begin(), conv_data.end(), conv.getData().begin());
	return conv;
}


OpenCLLocalConvolver::OpenCLLocalConvolver(std::shared_ptr<OpenCL_env> opencl_env) :
	env(opencl_env)
{
	if (!env) {
		throw invalid_parameter("Empty OpenCL environment given to OpenCLLocalConvolver");
	}
}

Image OpenCLLocalConvolver::convolve(const Image &src, const Image &krn, const Mask &mask)
{
	try {
		return _convolve(src, krn, mask);
	} catch (const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error while convolving: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

Image OpenCLLocalConvolver::_convolve(const Image &src, const Image &krn, const Mask &mask) {

	// We use a group size of 16x16, so let's extend the src image
	// to the next multiple of 16
	auto clpad_x = (16 - (src.getWidth() % 16)) % 16;
	auto clpad_y = (16 - (src.getHeight() % 16)) % 16;
	const Image clpad_src = src.extend(src.getWidth() + clpad_x,
	                                   src.getHeight() + clpad_y, 0, 0);

	// Convolve using the appropriate data type
	Image result;
	if (env->use_double) {
		result = _clpadded_convolve<double>(clpad_src, krn, src);
	}
	else {
		result = _clpadded_convolve<float>(clpad_src, krn, src);
	}

	// Crop the resulting image, mask
	return result.crop(src.getWidth(), src.getHeight(), 0, 0) & mask;
}

template<typename T>
Image OpenCLLocalConvolver::_clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src) {

	using cl::Buffer;
	using cl::Event;
	using cl::Kernel;
	using cl::Local;
	using cl::NDRange;
	using cl::NullRange;

	Buffer src_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, src.getSize());
	Buffer krn_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, krn.getSize());
	Buffer conv_buf = env->get_buffer<T>(CL_MEM_WRITE_ONLY, src.getSize());

	std::vector<T> src_data(src.getSize());
	std::copy(src.getData().begin(), src.getData().end(), src_data.begin());
	std::vector<T> krn_data(krn.getSize());
	std::copy(krn.getData().begin(), krn.getData().end(), krn_data.begin());

	// Write both images' data to the device
	Event src_wevt = env->queue_write(src_buf, src_data.data());
	Event krn_wevt = env->queue_write(krn_buf, krn_data.data());

	// We need this much local memory on each local group
	auto local_size = sizeof(T);
	local_size *= (16 + 2 * (krn.getWidth() / 2));
	local_size *= (16 + 2 * (krn.getHeight() / 2));

	if (env->max_local_memory() < local_size) {
		std::ostringstream os;
		os << "Not enough local memory available for OpenCL local 2D convolution. ";
		os << "Required: " << local_size << ", available: " << env->max_local_memory();
		throw opencl_error(os.str());
	}

	// Prepare the kernel
	auto kname = std::string("convolve_local_") + float_traits<T>::name;
	Kernel clKernel = env->get_kernel(kname);
	clKernel.setArg(0, src_buf);
	clKernel.setArg(1, orig_src.getWidth());
	clKernel.setArg(2, orig_src.getHeight());
	clKernel.setArg(3, krn_buf);
	clKernel.setArg(4, krn.getWidth());
	clKernel.setArg(5, krn.getHeight());
	clKernel.setArg(6, conv_buf);
	clKernel.setArg(7, Local(local_size));

	// Execute
	std::vector<Event> exec_wait_evts {src_wevt, krn_wevt};
	auto exec_evt = env->queue_kernel(clKernel, NDRange(src.getWidth(), src.getHeight()), &exec_wait_evts, NDRange(16, 16));

	// Read and good bye
	std::vector<T> conv_data(src.getSize());
	std::vector<Event> read_wait_evts {exec_evt};
	Event read_evt = env->queue_read(conv_buf, conv_data.data(), &read_wait_evts);
	read_evt.wait();

	Image conv(src.getWidth(), src.getHeight());
	std::copy(conv_data.begin(), conv_data.end(), conv.getData().begin());
	return conv;
}


#endif // PROFIT_OPENCL

std::shared_ptr<Convolver> create_convolver(const ConvolverType type, const ConvolverCreationPreferences &prefs)
{
	switch(type) {
		case BRUTE:
			return std::make_shared<BruteForceConvolver>();
#ifdef PROFIT_OPENCL
		case OPENCL:
			return std::make_shared<OpenCLConvolver>(prefs.opencl_env);
		case OPENCL_LOCAL:
			return std::make_shared<OpenCLLocalConvolver>(prefs.opencl_env);
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
		case FFT:
			return std::make_shared<FFTConvolver>(prefs.src_width, prefs.src_height,
			                                      prefs.krn_width, prefs.krn_height,
			                                      prefs.effort, prefs.plan_omp_threads,
			                                      prefs.reuse_krn_fft);
#endif // PROFIT_FFTW
		default:
			// Shouldn't happen
			throw invalid_parameter("Unsupported convolver type: " + std::to_string(type));
	}
}

std::shared_ptr<Convolver> create_convolver(const std::string &type, const ConvolverCreationPreferences &prefs)
{
	if (type == "brute") {
		return create_convolver(BRUTE, prefs);
	}
#ifdef PROFIT_OPENCL
	else if (type == "opencl") {
		return create_convolver(OPENCL, prefs);
	}
	else if (type == "opencl-local") {
		return create_convolver(OPENCL_LOCAL, prefs);
	}
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
	else if (type == "fft") {
		return create_convolver(FFT, prefs);
	}
#endif // PROFIT_FFTW

	std::ostringstream os;
	os << "Convolver of type " << type << " is not supported";
	throw invalid_parameter(os.str());
}


} /* namespace profit */
