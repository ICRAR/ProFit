/**
 * Public header file for the image convolution classes and methods
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
#include "profit/fft.h"
#include "profit/image.h"
#include "profit/opencl.h"

namespace profit
{

/**
 * The types of convolvers supported by libprofit
 */
enum ConvolverType {

	/// @copydoc BruteForceConvolver
	BRUTE_OLD = 0,

	/// @copydoc AssociativeBruteForceConvolver
	BRUTE,

	/// @copydoc OpenCLConvolver
	OPENCL,

	/// @copydoc FFTConvolver
	FFT,

	/// @copydoc NULLConvolver
	NO_OP
};

/**
 * A convolver object convolves two images.
 *
 * This is the base class for all Convolvers. Deriving classes must implement
 * the convolve method, which performs the actual operation.
 */
class PROFIT_API Convolver {

public:
	static Point NO_OFFSET;

	virtual ~Convolver();

	/**
	 * Convolves image `src` with the kernel `krn`.
	 * A mask parameter also controls which pixels from the original image
	 * should be convolved. If empty, all pixels are convolved.
	 *
	 * If the convolver extends the original image to perform the convolution,
	 * users might want to have the extended image returned, instead of getting
	 * a cropped image (that will be the same size as `src`). This behaviour is
	 * controlled with the `crop` parameter. If the image is not cropped, the
	 * offset of the otherwise cropped result with respect to the uncropped one
	 * is optionally stored in offset_out.
	 *
	 * @param src The source image
	 * @param krn The convolution kernel
	 * @param mask An mask indicating which pixels of the resulting image should
	 *             be convolved
	 * @param crop If ``true`` return an image with the same dimensions of ``src``.
	 *             If ``false`` the image returned might be potentially bigger,
	 *             depending on the internal workings of the convolver.
	 * @param offset_out If `crop` is ``false`` and ``offset`` is different from
	 *             NO_OFFSET, stores the potential offset of the original image
	 *             with respect to the uncropped image returned by this method.
	 * @return The convolved image, optionally without the cropping caused due
	 *         to internal implementation details of the convolver. The potential
	 *         offset is written into offset_out.
	 */
	Image convolve(const Image &src, const Image &krn, const Mask &mask,
	               bool crop = true, Point &offset_out = NO_OFFSET);

	/**
	 * Returns the amount of padding that would be introduced by this convolver
	 * when convolving an image and a kernel of sizes @p src_dims and @p
	 * krn_dims, respectively. The padding is returned as a pair of points (or
	 * dimensions) representing the padding in two dimensions at the bottom and
	 * top ends of the result, respectively.
	 *
	 * This method does @b not perform any convolution; it simply calculates
	 * what @em would be the padding of the image with respect to the
	 * convolution result.
	 *
	 * @param src_dims The dimensions of the source image
	 * @param krn_dims The dimensions of the kernel
	 * @return The padding that the convolution process would introduce to the
	 * returned image
	 */
	virtual PointPair padding(const Dimensions &src_dims, const Dimensions &krn_dims) const;

protected:
	// Implemented by subclasses and called by convolve
	virtual
	Image convolve_impl(const Image &src, const Image &krn, const Mask &mask,
	                    bool crop = true, Point &offset_out = NO_OFFSET) = 0;

	Image mask_and_crop(Image &img, const Mask &mask, bool crop,
	                    const Dimensions &orig_dims, const Dimensions &ext_dims,
	                    const Point &ext_offset, Point &offset_out);
};

///
/// A set of preferences used to create convolvers.
///
class PROFIT_API ConvolverCreationPreferences {

public:
	ConvolverCreationPreferences() :
		src_dims(),
		krn_dims(),
		omp_threads(1),
		opencl_env(),
		effort(effort_t::ESTIMATE),
		reuse_krn_fft(false),
		instruction_set(simd_instruction_set::AUTO)
	{};

	ConvolverCreationPreferences(
	    Dimensions src_dims, Dimensions krn_dims, unsigned int omp_threads,
	    OpenCLEnvPtr opencl_env, effort_t effort, bool reuse_krn_fft,
	    simd_instruction_set instruction_set) :
		src_dims(src_dims),
		krn_dims(krn_dims),
		omp_threads(omp_threads),
		opencl_env(opencl_env),
		effort(effort),
		reuse_krn_fft(reuse_krn_fft),
		instruction_set(instruction_set)
	{};


	/// The dimensions of the image being convolved.
	Dimensions src_dims;

	/// The dimensions of the convolution kernel.
	Dimensions krn_dims;

	/// The amount of OpenMP threads (if OpenMP is available) to use by the
	/// convolver. Used by the FFT convolver (to create and execute the plan
	/// using OpenMP, when available) and the brute-force convolvers.
	unsigned int omp_threads;

	/// A pointer to an OpenCL environment. Used by the OPENCL convolvers.
	OpenCLEnvPtr opencl_env;

	/// The amount of effort to put into the plan creation. Used by the @ref FFT convolver.
	effort_t effort;

	/// Whether to reuse or not the FFT'd kernel or not. Used by the @ref FFT convolver.
	bool reuse_krn_fft;

	/// The extended instruction set to use. Used by the @ref BRUTE convolver
	simd_instruction_set instruction_set;
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
PROFIT_API ConvolverPtr
create_convolver(const ConvolverType type,
                 const ConvolverCreationPreferences &prefs = ConvolverCreationPreferences());

/**
 * Like create_convolver(ConvolverType, const ConvolverCreationPreferences &),
 * but indicating the convolver type as a string.
 *
 * @overload
 */
PROFIT_API ConvolverPtr
create_convolver(const std::string &type,
                 const ConvolverCreationPreferences &prefs = ConvolverCreationPreferences());

} /* namespace profit */

#endif /* PROFIT_CONVOLVE_H */