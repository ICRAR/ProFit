/**
 * Private header file for the image convolution classes implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
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

#include "profit/convolve.h"
#include "profit/fft_impl.h"
#include "profit/opencl_impl.h"

namespace profit {

/**
 * A brute-force convolver. It optionally uses OpenMP to accelerate the
 * convolution.
 */
class BruteForceConvolver : public Convolver {

public:
	explicit BruteForceConvolver(unsigned int omp_threads) :
		omp_threads(omp_threads) {}

	Image convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop = true, Point &offset_out = NO_OFFSET) override;

private:
	unsigned int omp_threads;
};

/**
 * A faster brute-force convolver. It optionally uses OpenMP to accelerate the
 * convolution.
 *
 * The difference between this and the BruteForceConvolver is that this
 * convolver explicitly states that the sums of the dot products that
 * make up the result of a single pixel are associative, and can be computed
 * separately, which enables better pipelining in most CPUs and thus faster
 * compute times (we have seen up to ~3x speedups). The result is not guaranteed
 * to be the exact same as the one coming from BruteForceConvolver. This is not
 * because one of them is mathematically incorrect (neither is actually), but
 * because IEEE floating-point math is not associative, and therefore different
 * operation sequences *might* yield different results.
 *
 * The internal loop structure of this class is also slightly different from
 * BruteForceConvolver, but is still pure CPU-based code.
 *
 * Additionally, and depending on the underlying CPU support, this convolver
 * can use dot product implementations based on SIMD operations available in
 * different CPU extended instruction sets. The default is to use the fastest
 * one available, although users might want to use a different one.
 */
template <simd_instruction_set SIMD>
class AssociativeBruteForceConvolver : public Convolver {

public:
	explicit AssociativeBruteForceConvolver(unsigned int omp_threads) : omp_threads(omp_threads) {};

protected:
	Image convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop = true, Point &offset_out = NO_OFFSET) override;

private:
	unsigned int omp_threads;
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
 * After convolution the result is cropped back (if required) to the original image's
 * dimensions starting at the center of the original image's mapping on the
 * extended image (i.e., ``(src_width/2, src_height/2)`` minus one if the
 * original dimensions are odd).
 *
 * This convolver has been implemented in such a way that no memory allocation
 * happens during convolution (other than the final Image's allocation) to
 * improve performance.
 */
class FFTConvolver : public Convolver {

public:
	explicit FFTConvolver(const Dimensions &src_dims, const Dimensions &krn_dims,
	             effort_t effort, unsigned int plan_omp_threads,
	             bool reuse_krn_fft);

	PointPair padding(const Dimensions &src_dims, const Dimensions &krn_dims) const override;

protected:
	Image convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop = true, Point &offset_out = NO_OFFSET) override;

private:

	void resize(const Dimensions &src_dims, const Dimensions &krn_dims);

	Point offset_after_convolution(const Dimensions &src_dims, const Dimensions &krn_dims) const;

	std::unique_ptr<FFTRealTransformer> fft_transformer;

	std::vector<std::complex<double>> src_fft;
	std::vector<std::complex<double>> krn_fft;
	Image ext_src;
	Image ext_krn;

	bool reuse_krn_fft;
	bool krn_fft_initialized;
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
	explicit OpenCLConvolver(OpenCLEnvImplPtr opencl_env);

	PointPair padding(const Dimensions &src_dims, const Dimensions &krn_dims) const override;

protected:
	Image convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop = true, Point &offset_out = NO_OFFSET) override;

private:
	OpenCLEnvImplPtr env;

	// returns the extra OpenCL-imposed padding
	Dimensions cl_padding(const Dimensions &src_dims) const;

	Image _convolve(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out);

	template<typename T>
	Image _clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src);
};

/**
 * Like OpenCLConvolver, but uses a local memory cache
 */
class OpenCLLocalConvolver : public Convolver {

public:
	explicit OpenCLLocalConvolver(OpenCLEnvImplPtr opencl_env);

protected:
	Image convolve_impl(const Image &src, const Image &krn, const Mask &mask, bool crop = true, Point &offset_out = NO_OFFSET) override;

private:
	OpenCLEnvImplPtr env;

	Image _convolve(const Image &src, const Image &krn, const Mask &mask, bool crop, Point &offset_out);

	template<typename T>
	Image _clpadded_convolve(const Image &src, const Image &krn, const Image &orig_src);
};

#endif // PROFIT_OPENCL

}