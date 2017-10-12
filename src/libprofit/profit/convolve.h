/**
 * Header file for the image convolution implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
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
#ifndef PROFIT_CONVOLVE_H
#define PROFIT_CONVOLVE_H

#include <complex>
#include <memory>
#include <vector>

#include "profit/common.h"
#include "profit/image.h"
#include "profit/fft.h"
#include "profit/opencl.h"

namespace profit
{

/**
 * The types of convolvers supported by libprofit
 */
enum ConvolverType {
	BRUTE = 0,
#ifdef PROFIT_OPENCL
	OPENCL,
	OPENCL_LOCAL,
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
	FFT,
#endif // PROFIT_FFTW
};


/**
 * A convolver object convolves two images.
 *
 * This is the base class for all Convolvers. Deriving classes must implement
 * the convolve method, which performs the actual operation.
 */
class Convolver {

public:
	virtual ~Convolver();

	/**
	 * Convolves image `src` with the kernel `krn`.
	 * A mask parameter also controls which pixels from the original image
	 * should be convolved. If empty, all pixels are convolved.

	 * @param src The source image
	 * @param krn The convolution kernel
	 * @param mask An mask indicating which pixels of the resulting image should
	 *             be convolved
	 * @return The convolved image
	 */
	virtual
	Image convolve(const Image &src, const Image &krn, const Mask &mask) = 0;

};

/**
 * A brute-force convolver.
 */
class BruteForceConvolver : public Convolver {

public:
	Image convolve(const Image &src, const Image &krn, const Mask &mask) override;

};

#ifdef PROFIT_FFTW

/**
 * A convolver that uses an FFTPlan to carry out FFT-based convolution.
 *
 * The result of the convolution of images im1 and im2 is::
 *
 *  res = iFFT(FFT(im1) * FFT(im2))
 *
 * To do this, this convolver creates extended versions of the input images.
 * The size of the new images is 4 times that of the source image, which is
 * assumed to be larger than the kernel. The extended version of the source
 * image contains the original image at (0,0), while the extended version of the
 * kernel image contains the original kernel centered at the original image's
 * new mapping (i.e., ``((src_width-krn_width)/2, (src_height-krn_height)/2)``).
 * After convolution the result is cropped back to the original image's
 * dimensions starting at the center of the original image's mapping on the
 * extended image (i.e., ``(src_width/2, src_height/2)`` minus one if the
 * original dimensions are odd).
 */
class FFTConvolver : public Convolver {

public:
	FFTConvolver(unsigned int src_width, unsigned int src_height,
	             unsigned int krn_width, unsigned int krn_height,
	             FFTPlan::effort_t effort, unsigned int plan_omp_threads,
	             bool reuse_krn_fft);

	Image convolve(const Image &src, const Image &krn, const Mask &mask) override;

private:
	std::unique_ptr<FFTPlan> plan;

	std::vector<std::complex<double>> krn_fft;

	bool reuse_krn_fft;
};

#endif /* PROFIT_FFTW */

#ifdef PROFIT_OPENCL

/**
 * A brute-force convolver that is implemented using OpenCL
 *
 * Depending on the floating-point support found at runtime in the given OpenCL
 * environment this convolver will use a float-based or a double-based kernel.
 */
class OpenCLConvolver : public Convolver {

public:
	OpenCLConvolver(OpenCLEnvPtr opencl_env);

	Image convolve(const Image &src, const Image &krn, const Mask &mask) override;

private:
	OpenCLEnvPtr env;

	Image _convolve(const Image &src, const Image &krn, const Mask &mask);

	template<typename T>
	Image _clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src);
};

/**
 * Like OpenCLConvolver, but uses a local memory cache
 */
class OpenCLLocalConvolver : public Convolver {

public:
	OpenCLLocalConvolver(OpenCLEnvPtr opencl_env);

	Image convolve(const Image &src, const Image &krn, const Mask &mask) override;

private:
	OpenCLEnvPtr env;

	Image _convolve(const Image &src, const Image &krn, const Mask &mask);

	template<typename T>
	Image _clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src);
};


#endif // PROFIT_OPENCL


///
/// A set of preferences used to create convolvers.
///
class ConvolverCreationPreferences {

public:
	ConvolverCreationPreferences() :
		src_width(0),
		src_height(0),
		krn_width(0),
		krn_height(0)
#ifdef PROFIT_OPENCL
		,opencl_env()
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
		,effort(FFTPlan::ESTIMATE)
		,plan_omp_threads()
		,reuse_krn_fft(false)
#endif // PROFIT_FFTW
	{};

	/// The width of the image being convolved.
	unsigned int src_width;

	/// The height of the image being convolved.
	unsigned int src_height;

	/// The width of the convolution kernel.
	unsigned int krn_width;

	/// The height of the convolution kernel.
	unsigned int krn_height;

#ifdef PROFIT_OPENCL
	/// A pointer to an OpenCL environment. Used by the OpenCL convolvers.
	OpenCLEnvPtr opencl_env;
#endif // PROFIT_OPENCL

#ifdef PROFIT_FFTW

	/// The amount of effort to put into the plan creation. Used by the FFT convolver.
	FFTPlan::effort_t effort;

	/// The amount of OpenMP threads (if OpenMP is available) to use to create
	/// and execute the plan. Used by the FFT convolver.
	unsigned int plan_omp_threads;

	/// Whether to reuse or not the FFT'd kernel or not. Used by the FFT convolver.
	bool reuse_krn_fft;
#endif // PROFIT_FFTW

};

/// Handy typedef for shared pointers to Convolver objects
typedef std::shared_ptr<Convolver> ConvolverPtr;

/**
 * Creates a new convolver of type `type` with preferences `prefs`
 *
 * @param type The type of convolver to create
 * @param prefs The creation preferences used to create the new convolver
 * @return A shared pointer to a new convolver
 */
ConvolverPtr
create_convolver(const ConvolverType type,
                 const ConvolverCreationPreferences &prefs = ConvolverCreationPreferences());

/**
 * Like create_convolver(ConvolverType, const ConvolverCreationPreferences &),
 * but indicating the convolver type as a string.
 *
 * @overload
 */
ConvolverPtr
create_convolver(const std::string &type,
                 const ConvolverCreationPreferences &prefs = ConvolverCreationPreferences());

} /* namespace profit */

#endif /* PROFIT_CONVOLVE_H */