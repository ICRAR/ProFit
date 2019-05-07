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

#include "profit/common.h"
#include "profit/convolver_impl.h"
#include "profit/exceptions.h"
#include "profit/utils.h"


namespace profit
{

Point Convolver::NO_OFFSET;

Convolver::~Convolver()
{
	// no-op
}

Image Convolver::convolve(const Image &src, const Image &krn, const Mask &mask,
                          bool crop, Point &offset_out)
{
	// convolve_impl requires images to be not smaller than kernels. If that's
	// not the case we expand the image, convolve, then crop down the expansion.
	// During the final cropping we make sure we only crop out the padding we
	// added, and not the extra padding the convolution implementation might
	// have added
	auto src_dims = src.getDimensions();
	auto krn_dims = krn.getDimensions();
	if (krn_dims.x > src_dims.x || krn_dims.y > src_dims.y) {
		auto new_dims = max(src_dims, krn_dims);
		auto dims_diff = new_dims - src_dims;
		auto extended_src = src.extend(new_dims);
		auto extended_mask = Mask{};
		if (mask) {
			extended_mask = mask.extend(new_dims);
		}
		auto result = convolve_impl(extended_src, krn, extended_mask, crop, offset_out);
		return result.crop(result.getDimensions() - dims_diff);
	}
	else {
		return convolve_impl(src, krn, mask, crop, offset_out);
	}
}

Image Convolver::mask_and_crop(Image &img, const Mask &mask, bool crop, const Dimensions orig_dims, const Dimensions &ext_dims, const Point &ext_offset, Point &offset_out) {

	// No cropping requested
	// Extend the mask to the size of the image, apply it,
	// and save the offset of the original image with respect to the extension
	// into offset_out
	if (!crop) {
		if (&offset_out != &NO_OFFSET) {
			offset_out = ext_offset;
		}
		if (mask) {
			img &= mask.extend(ext_dims, ext_offset);
		}
		return img;
	}

	// Return cropped, after applying the mask
	return img.crop(orig_dims, ext_offset) & mask;
}

Image BruteForceConvolver::convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out)
{
	const auto src_dims = src.getDimensions();
	const auto krn_dims = krn.getDimensions();
	const auto src_width = src_dims.x;
	const auto src_height = src_dims.y;
	const auto krn_width = krn_dims.x;
	const auto krn_height = krn_dims.y;

	const unsigned int krn_half_width = krn_width / 2;
	const unsigned int krn_half_height = krn_height / 2;

	Image convolution(src_dims);

	auto krn_end = krn.cend();

	/* Convolve! */
	/* Loop around the output image first... */
#ifndef PROFIT_OPENMP
	UNUSED(omp_threads);
#else
	bool use_omp = omp_threads > 1;
	#pragma omp parallel for collapse(2) schedule(dynamic, 10) if(use_omp) num_threads(omp_threads)
#endif // PROFIT_OPENMP
	for (unsigned int j = 0; j < src_height; j++) {
		for (unsigned int i = 0; i < src_width; i++) {

			auto im_idx = i + j * src_width;

			/* Don't convolve this pixel */
			if( mask && !mask[im_idx]) {
				convolution[im_idx] = 0;
				continue;
			}

			double pixel = 0;
			auto krnPtr = krn_end - 1;
			auto srcPtr2 = src.begin() + im_idx - krn_half_width - krn_half_height*src_width;

			/* ... now loop around the kernel */
			for (unsigned int l = 0; l < krn_height; l++) {

				int src_j = (int)j + (int)l - (int)krn_half_height;
				for (unsigned int k = 0; k < krn_width; k++) {

					int src_i = (int)i + (int)k - (int)krn_half_width;

					if( src_i >= 0 && (unsigned int)src_i < src_width &&
					    src_j >= 0 && (unsigned int)src_j < src_height ) {
						pixel +=  *srcPtr2 * *krnPtr;
					}

					srcPtr2++;
					krnPtr--;
				}
				srcPtr2 += src_width - krn_width;
			}

			convolution[im_idx] = pixel;
		}
	}

	return convolution;
}


Image AssociativeBruteForceConvolver::convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out)
{

	const auto src_dims = src.getDimensions();
	const auto krn_dims = krn.getDimensions();
	const auto src_width = src_dims.x;
	const auto src_height = src_dims.y;
	const auto krn_width = krn_dims.x;
	const auto krn_height = krn_dims.y;

	const unsigned int krn_half_width = krn_width / 2;
	const unsigned int krn_half_height = krn_height / 2;

	Image convolution(src_dims);

	const size_t src_krn_offset = krn_half_width + krn_half_height*src_width;
	const auto src_skip = src_width - krn_width;

	/* Convolve! */
	/* Loop around the output image first... */
#ifndef PROFIT_OPENMP
	UNUSED(omp_threads);
#else
	bool use_omp = omp_threads > 1;
	#pragma omp parallel for collapse(2) schedule(dynamic, 10) if(use_omp) num_threads(omp_threads)
#endif // PROFIT_OPENMP
	for (unsigned int j = 0; j < src_height; j++) {
		for (unsigned int i = 0; i < src_width; i++) {

			auto im_idx = i + j * src_width;

			/* Don't convolve this pixel */
			if (mask && !mask[im_idx]) {
				convolution[im_idx] = 0;
				continue;
			}

			double pixel = 0;

			size_t krnPtr = krn.size() - 1;
			size_t srcPtr2 = im_idx;
			bool suboffset = false;

			unsigned int l_min = 0;
			unsigned int l_max = krn_height;
			unsigned int l_incr = 0;

			if (j < krn_half_height) {
				l_min = krn_half_height - j;
				srcPtr2 += l_min * src_width;
				krnPtr -= l_min * krn_width;
			}
			else if ((j + krn_half_height) >= src_height) {
				// TODO: maybe shouldn't be an else if we support krn > img size?
				l_max = src_height + krn_half_height - j;
				l_incr = krn_height - l_max;
			}

			for (size_t l = l_min; l < l_max; l++) {

				unsigned int k_min = 0;
				unsigned int k_max = krn_width;
				unsigned int k_incr = 0;

				if (i < krn_half_width) {
					k_min = krn_half_width - i;
					srcPtr2 += k_min;
					krnPtr -= k_min;
				}
				else if ((i + krn_half_width) >= src_width)
				{
					// TODO: maybe shouldn't be an else-if if we support krn > img size?
					k_max = src_width + krn_half_width - i;
					k_incr = krn_width - k_max;
				}

				if (!suboffset && srcPtr2 >= src_krn_offset)
				{
					srcPtr2 -= src_krn_offset;
					suboffset = true;
				}
				const size_t k_n = k_max - k_min;

				// Sum multiplications first, then add up to pixel.
				// This means we explicitly tell the compiler that:
				//
				//  a + b + c + d == (a + b + c) + d
				//
				// By default floating point arithmetic is not associative,
				// and therefore the compiler will not create the temporary
				// "buf" variable, unless compiling with -ffast-math et al.
				// Doing this buffering allows compilers to use an extra
				// register, which in turn yields better instruction pipelining.
				double buf = 0;

				// On top of the associativity described above,
				// we also manually unroll the for loop into four separate
				// multiply-add operations. allows compilers to optimize even
				// further, because there is more explicit associativity and
				// thus better pipelining
				//
				// Also, note that clang needs an explicit -ffp-contract=fast
				// to generate fused multiply-add instructions (which gcc does
				// for default). This is not only important here, but also in
				// the original version of our convolution method.
				//
				// TODO: The generated SSE/AVX instructions are still not
				//       vectorized (e.g., vfmaddsd instead of vfmaddpd). This
				//       is because the compiler cannot guarantee the alignment
				//       of the arrays. The difficulty on doing that lies on the
				//       the fact that both arrays move separately, so it's
				//       difficult to make that bring that kind of assurance
				//       (other than copying data to an aligned buffer).
				//       An additional benefit from generating vectorized
				//       instructions is that the compiler can fully use the
				//       YMM registers (and not only half of the XMM registers,
				//       as it is doing now) leading to yet better performance.
				for (size_t k = 0; k < k_n / 4; k++) {
					double tmp1 = src[srcPtr2 + k * 4]     * krn[krnPtr - k * 4];
					double tmp2 = src[srcPtr2 + k * 4 + 1] * krn[krnPtr - k * 4 - 1];
					double tmp3 = src[srcPtr2 + k * 4 + 2] * krn[krnPtr - k * 4 - 2];
					double tmp4 = src[srcPtr2 + k * 4 + 3] * krn[krnPtr - k * 4 - 3];
					buf += (tmp1 + tmp3) + (tmp2 + tmp4);
				}
				switch (k_n % 4) {
					case 3:
						buf += src[srcPtr2 + k_n - 3] * krn[krnPtr - k_n + 3] + \
						       src[srcPtr2 + k_n - 2] * krn[krnPtr - k_n + 2] + \
						       src[srcPtr2 + k_n - 1] * krn[krnPtr - k_n + 1];
						break;

					case 2:
						buf += src[srcPtr2 + k_n - 2] * krn[krnPtr - k_n + 2] + \
						       src[srcPtr2 + k_n - 1] * krn[krnPtr - k_n + 1];
						break;

					case 1:
						buf += src[srcPtr2 + k_n - 1] * krn[krnPtr - k_n + 1];
						break;

					case 0:
						break;
				}

				pixel += buf;
				srcPtr2 += k_n;
				krnPtr -= k_n;

				srcPtr2 += k_incr;
				krnPtr -= k_incr;
				srcPtr2 += src_skip;
			}

			srcPtr2 += l_incr * krn_width;
			krnPtr -= l_incr * krn_width;

			convolution[im_idx] = pixel;
		}
	}

	return convolution;

}

#ifdef PROFIT_FFTW
FFTConvolver::FFTConvolver(const Dimensions &src_dims, const Dimensions &krn_dims,
                           effort_t effort, unsigned int plan_omp_threads,
                           bool reuse_krn_fft) :
	fft_transformer(),
	krn_fft(),
	reuse_krn_fft(reuse_krn_fft)
{
	auto effective_dims = max(src_dims, krn_dims);
	auto convolution_size = 4 * effective_dims.x * effective_dims.y;
	fft_transformer = std::unique_ptr<FFTTransformer>(new FFTRealTransformer(convolution_size, effort, plan_omp_threads));
}

Image FFTConvolver::convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out)
{

	typedef std::complex<double> complex;

	auto src_dims = src.getDimensions();
	auto krn_dims = krn.getDimensions();

	// Create extended images first
	auto ext_dims = src_dims * 2;
	Image ext_img = src.extend(ext_dims);

	// Forward FFTs
	std::vector<complex> src_fft = fft_transformer->forward(ext_img);
	if (krn_fft.empty()) {
		auto krn_start = (src_dims - krn_dims) / 2;
		Image ext_krn = krn.extend(ext_dims, krn_start);
		krn_fft = fft_transformer->forward(ext_krn);
	}

	// element-wise multiplication
	std::transform(src_fft.begin(), src_fft.end(), krn_fft.begin(), src_fft.begin(), std::multiplies<complex>());

	if (!reuse_krn_fft) {
		krn_fft.clear();
	}

	// inverse FFT and scale down
	Image res(fft_transformer->backward(src_fft), ext_dims);
	res /= res.size();

	// The resulting image now starts at x_offset/y_offset
	// even image and odd kernel requires slight adjustment
	auto ext_offset = src_dims / 2;
	if (src_dims.x % 2 == 0 || krn_dims.x % 2 == 0) {
		ext_offset.x -= 1;
	}
	if (src_dims.y % 2 == 0 || krn_dims.y % 2 == 0) {
		ext_offset.y -= 1;
	}

	return mask_and_crop(res, mask, crop, src_dims, ext_dims, ext_offset, offset_out);
}

#endif /* PROFIT_FFTW */

#ifdef PROFIT_OPENCL
OpenCLConvolver::OpenCLConvolver(OpenCLEnvImplPtr opencl_env) :
	env(opencl_env)
{
	if (!env) {
		throw invalid_parameter("Empty OpenCL environment given to OpenCLConvolver");
	}
}

Image OpenCLConvolver::convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out)
{
	try {
		return _convolve(src, krn, mask, crop, offset_out);
	} catch (const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error while convolving: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

Image OpenCLConvolver::_convolve(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out) {

	// We use a group size of 16x16, so let's extend the src image
	// to the next multiple of 16
	auto src_dims = src.getDimensions();
	auto clpad_dims = (16 - (src_dims % 16)) % 16;
	auto ext_dims = src.getDimensions() + clpad_dims;
	const Image clpad_src = src.extend(ext_dims);

	// Convolve using the appropriate data type
	Image result;
	if (env->use_double) {
		result = _clpadded_convolve<double>(clpad_src, krn, src);
	}
	else {
		result = _clpadded_convolve<float>(clpad_src, krn, src);
	}

	return mask_and_crop(result, mask, crop, src_dims, ext_dims, {0, 0}, offset_out);
}

template<typename T>
Image OpenCLConvolver::_clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src) {

	using cl::Buffer;
	using cl::Event;
	using cl::Kernel;
	using cl::NDRange;
	using cl::NullRange;

	Buffer src_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, src.size());
	Buffer krn_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, krn.size());
	Buffer conv_buf = env->get_buffer<T>(CL_MEM_WRITE_ONLY, src.size());

	std::vector<T> src_data(src.size());
	std::copy(src.begin(), src.end(), src_data.begin());
	std::vector<T> krn_data(krn.size());
	std::copy(krn.begin(), krn.end(), krn_data.begin());

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
	auto exec_evt = env->queue_kernel(clKernel, NDRange(src.getWidth(), src.getHeight()), &exec_wait_evts);

	// Read and good bye
	std::vector<Event> read_wait_evts {exec_evt};
	std::vector<T> conv_data(src.size());
	Event read_evt = env->queue_read(conv_buf, conv_data.data(), &read_wait_evts);
	read_evt.wait();

	Image conv(src.getDimensions());
	std::copy(conv_data.begin(), conv_data.end(), conv.begin());
	return conv;
}


OpenCLLocalConvolver::OpenCLLocalConvolver(OpenCLEnvImplPtr opencl_env) :
	env(opencl_env)
{
	if (!env) {
		throw invalid_parameter("Empty OpenCL environment given to OpenCLLocalConvolver");
	}
}

Image OpenCLLocalConvolver::convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out)
{
	try {
		return _convolve(src, krn, mask, crop, offset_out);
	} catch (const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error while convolving: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

Image OpenCLLocalConvolver::_convolve(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out) {

	// We use a group size of 16x16, so let's extend the src image
	// to the next multiple of 16
	auto src_dims = src.getDimensions();
	auto clpad_dims = (16 - src_dims % 16) % 16;
	auto ext_dims = src_dims + clpad_dims;
	const Image clpad_src = src.extend(ext_dims);

	// Convolve using the appropriate data type
	Image result;
	if (env->use_double) {
		result = _clpadded_convolve<double>(clpad_src, krn, src);
	}
	else {
		result = _clpadded_convolve<float>(clpad_src, krn, src);
	}

	return mask_and_crop(result, mask, crop, src_dims, ext_dims, {0, 0}, offset_out);
}

template<typename T>
Image OpenCLLocalConvolver::_clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src) {

	using cl::Buffer;
	using cl::Event;
	using cl::Kernel;
	using cl::Local;
	using cl::NDRange;
	using cl::NullRange;

	Buffer src_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, src.size());
	Buffer krn_buf = env->get_buffer<T>(CL_MEM_READ_ONLY, krn.size());
	Buffer conv_buf = env->get_buffer<T>(CL_MEM_WRITE_ONLY, src.size());

	std::vector<T> src_data(src.size());
	std::copy(src.begin(), src.end(), src_data.begin());
	std::vector<T> krn_data(krn.size());
	std::copy(krn.begin(), krn.end(), krn_data.begin());

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
	std::vector<T> conv_data(src.size());
	std::vector<Event> read_wait_evts {exec_evt};
	Event read_evt = env->queue_read(conv_buf, conv_data.data(), &read_wait_evts);
	read_evt.wait();

	Image conv(src.getDimensions());
	std::copy(conv_data.begin(), conv_data.end(), conv.begin());
	return conv;
}


#endif // PROFIT_OPENCL

ConvolverPtr create_convolver(const ConvolverType type, const ConvolverCreationPreferences &prefs)
{
	switch(type) {
		case BRUTE_OLD:
			return std::make_shared<BruteForceConvolver>(prefs.omp_threads);
		case BRUTE:
			return std::make_shared<AssociativeBruteForceConvolver>(prefs.omp_threads);
#ifdef PROFIT_OPENCL
		case OPENCL:
			return std::make_shared<OpenCLConvolver>(OpenCLEnvImpl::fromOpenCLEnvPtr(prefs.opencl_env));
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
		case FFT:
			return std::make_shared<FFTConvolver>(prefs.src_dims, prefs.krn_dims,
			                                      prefs.effort, prefs.omp_threads,
			                                      prefs.reuse_krn_fft);
#endif // PROFIT_FFTW
		default:
			// Shouldn't happen
			throw invalid_parameter("Unsupported convolver type: " + std::to_string(type));
	}
}

ConvolverPtr create_convolver(const std::string &type, const ConvolverCreationPreferences &prefs)
{
	if (type == "brute-old") {
		return create_convolver(BRUTE_OLD, prefs);
	}
	else if (type == "brute") {
		return create_convolver(BRUTE, prefs);
	}
#ifdef PROFIT_OPENCL
	else if (type == "opencl") {
		return create_convolver(OPENCL, prefs);
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
