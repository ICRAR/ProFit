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
	const auto src_width = src.getWidth();
	const auto src_height = src.getHeight();
	const auto krn_width = krn.getWidth();
	const auto krn_height = krn.getHeight();

	const unsigned int krn_half_width = krn_width / 2;
	const unsigned int krn_half_height = krn_height / 2;

	Image convolution(src_width, src_height);

	const auto &src_data = src.getData();
	auto &out = convolution.getData();
	const auto &mask_data = mask.getData();
	const auto &krn_data = krn.getData();
	const size_t src_krn_offset = krn_half_width + krn_half_height*src_width;
	const auto src_skip = src_width - krn_width;

	/* Convolve! */
	/* Loop around the output image first... */
#ifdef PROFIT_OPENMP
	bool use_omp = omp_threads > 1;
	#pragma omp parallel for collapse(2) schedule(dynamic, 10) if(use_omp) num_threads(omp_threads)
#endif // PROFIT_OPENMP
	for (unsigned int j = 0; j < src_height; j++) {
		for (unsigned int i = 0; i < src_width; i++) {

			auto im_idx = i + j * src_width;

			/* Don't convolve this pixel */
			if( !mask.empty() and !mask_data[im_idx]) {
				out[im_idx] = 0;
				continue;
			}

			double pixel = 0;

			size_t krnPtr = krn_data.size() - 1;
			size_t srcPtr = im_idx;
			bool suboffset = false;

			unsigned int l_min = 0;
			unsigned int l_max = krn_height;
			unsigned int l_incr = 0;

			if(j < krn_half_height) {
				l_min = krn_half_height - j;
				srcPtr += l_min * src_width;
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
					srcPtr += k_min;
					krnPtr -= k_min;
				}
				else if((i + krn_half_width) >= src_width)
				{
					// TODO: maybe shouldn't be an else-if if we support krn > img size?
					k_max = src_width + krn_half_width - i;
					k_incr = krn_width - k_max;
				}

				if (!suboffset and srcPtr >= src_krn_offset)
				{
					srcPtr -= src_krn_offset;
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
				
				const auto K_MAX = k_n / 4;
				for (size_t k = 0; k < K_MAX; k++) {
					double tmp1 = src_data[srcPtr + k*4] * krn_data[krnPtr - k*4];
					double tmp2 = src_data[srcPtr + k*4 + 1] * krn_data[krnPtr - k*4 - 1];
					double tmp3 = src_data[srcPtr + k*4 + 2] * krn_data[krnPtr - k*4 - 2];
					double tmp4 = src_data[srcPtr + k*4 + 3] * krn_data[krnPtr - k*4 - 3];
					buf += (tmp1 + tmp3) + (tmp2 + tmp4);
				}
				switch (k_n % 4) {
					case 3:
						buf += src_data[srcPtr + k_n - 3] * krn_data[krnPtr - k_n + 3] + \
						       src_data[srcPtr + k_n - 2] * krn_data[krnPtr - k_n + 2] + \
						       src_data[srcPtr + k_n - 1] * krn_data[krnPtr - k_n + 1];
						break;

					case 2:
						buf += src_data[srcPtr + k_n - 2] * krn_data[krnPtr - k_n + 2] + \
						       src_data[srcPtr + k_n - 1] * krn_data[krnPtr - k_n + 1];
						break;

					case 1:
						buf += src_data[srcPtr + k_n - 1] * krn_data[krnPtr - k_n + 1];
						break;

					case 0:
						break;
				}

				pixel += buf;
				srcPtr += k_n;
				krnPtr -= k_n;

				srcPtr += k_incr;
				krnPtr -= k_incr;
				srcPtr += src_skip;
			}

			srcPtr += l_incr * krn_width;
			krnPtr -= l_incr * krn_width;

			out[im_idx] = pixel;
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

Image FFTConvolver::convolve(const Image &src, const Image &krn, const Mask &mask, 
  const Image * extra, bool extraissrc)
{
  const bool hasextra = extra != nullptr;
  
	typedef std::complex<double> complex;

	const auto src_width = src.getWidth();
	const auto src_height = src.getHeight();
	const auto krn_width = krn.getWidth();
	const auto krn_height = krn.getHeight();
	const auto extra_width = hasextra ? extra->getWidth() : 0;
	const auto extra_height = hasextra ? extra->getHeight() : 0;
	const bool extra_src_dim_same = !hasextra ? false :
	  (extra_width == src_width && extra_height == src_height);

	if(hasextra)
	{
	  if(!(extra_width <= src_width) && !(extra_height <= src_height))
	  {
	    // TODO: Output more info like actual dimensions
	    throw fft_error("Error! Extra image/kernel larger than source image.");
	  }
	}
	
	// Zero padding
	auto ext_width = 2 * src_width;
	auto ext_height = 2 * src_height;

	// Forward FFTs
	std::vector<complex> src_fft;
	if(hasextra)
	{
	  src_fft.resize(ext_width*ext_height);

	  const auto & srcvec = src.getData();
	  const auto & extravec = extra->getData();

	  // TODO: Test to see if this is actually faster.
	  if(extra_src_dim_same)
	  {
	    for(unsigned int j = 0; j < src_height; j++) {
  		  for(unsigned int i = 0; i < src_width; i++) {
  			  src_fft[i + j*ext_width].real(srcvec[i + j*src_width]);
  		    src_fft[i + j*ext_width].imag(extravec[i + j*src_width]);
  		  }
  	  }
	  } else {
	    auto start_x = 0;
	    auto start_y = 0;
  	  for(unsigned int j = 0; j < src_height; j++) {
  		  for(unsigned int i = 0; i < src_width; i++) {
  			  src_fft[(i+start_x) + (j+start_y)*ext_width].real(srcvec[i + j*src_width]);
  		  }
  	  }
  	  // Shift the image a bit to allow for some zero padding which won't be discarded
  	  start_x = (src_width-extra_width)/2;
  	  start_y = (src_height-extra_height)/2;
  	  for(unsigned int j = 0; j < extra_height; j++) {
  		  for(unsigned int i = 0; i < extra_width; i++) {
  			  src_fft[(i+start_x) + (j+start_y)*ext_width].imag(extravec[i + j*extra_width]);
  		  }
  	  }
	  }
	  // TODO: Can this be done in-place?
	  src_fft = plan->forward(src_fft);
	} else {
	  // Create extended images first
	  Image ext_img = src.extend(ext_width, ext_height, 0, 0);
	  src_fft = plan->forward(ext_img);
	}
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

	// inverse FFT and scale down
	Image cropped;
	if(hasextra)
	{
	  std::vector<complex> res = plan->backward(src_fft);
	  const auto src_height2 = 2*src_height;
	  std::vector<double> cropvec(src_width*src_height2);
	  
    // The result is an image of twice the original width, for (some) convenience
    // It keeps the function signature the same, at least, although it then
    // requires subsetting on the user side
    
    for(unsigned int j = 0; j < src_height; j++) {
		  for(unsigned int i = 0; i < src_width; i++) {
			  cropvec[i + j*src_width] = res[i + x_offset + (j + y_offset)*ext_width].real();
		    cropvec[i + j*src_width + src_width*src_height] = res[i + x_offset + (j + y_offset)*ext_width].imag();
		  }
	  }
    if(!mask.empty())
    {
      const auto & maskvec = mask.getData();
      std::transform(cropvec.begin(), cropvec.begin() + src_width*src_height,
        maskvec.begin(), cropvec.begin(),
    		[](const double i, const bool m) {
    			return m ? i : 0.;
    	});
      std::transform(cropvec.begin() + src_width*src_height, cropvec.end(),
        maskvec.begin(), cropvec.begin() + src_width*src_height,
    		[](const double i, const bool m) {
    			return m ? i : 0.;
    	});
    }
    cropped = Image(cropvec, src_width, src_height2);
	} else {
	  Image res = Image(plan->backward_real(src_fft), ext_width, ext_height);
	  cropped = res.crop(src_width, src_height, x_offset, y_offset);
	  cropped &= mask;
	}
	cropped /= ext_width*ext_height;
	return cropped;
}

#endif /* PROFIT_FFTW */

#ifdef PROFIT_OPENCL
OpenCLConvolver::OpenCLConvolver(OpenCLEnvPtr opencl_env) :
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
	auto exec_evt = env->queue_kernel(clKernel, NDRange(src.getWidth(), src.getHeight()), &exec_wait_evts);

	// Read and good bye
	std::vector<Event> read_wait_evts {exec_evt};
	std::vector<T> conv_data(src.getSize());
	Event read_evt = env->queue_read(conv_buf, conv_data.data(), &read_wait_evts);
	read_evt.wait();

	Image conv(src.getWidth(), src.getHeight());
	std::copy(conv_data.begin(), conv_data.end(), conv.getData().begin());
	return conv;
}


OpenCLLocalConvolver::OpenCLLocalConvolver(OpenCLEnvPtr opencl_env) :
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

ConvolverPtr create_convolver(const ConvolverType type, const ConvolverCreationPreferences &prefs)
{
	switch(type) {
		case BRUTE:
			return std::make_shared<BruteForceConvolver>(prefs.omp_threads);
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
