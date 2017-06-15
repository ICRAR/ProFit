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
#include "profit/fft.h"

namespace profit
{

/**
 * Brute-force convolves image `src` with the kernel `krn`.
 *
 * A mask parameter also controls which pixels from the original image should be
 * convolved. If NULL all pixels are convolved.

 * @param src The source image
 * @param src_width The source image's width
 * @param src_height The source image's height
 * @param krn The convolution kernel
 * @param krn_width The kernel's width
 * @param krn_height The kernel's height
 * @param mask An optional boolean mask indicating which pixels of the resulting
 *             image should be convolved
 * @return The convolved image
 */
std::vector<double>
convolve(const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
         const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
         const std::vector<bool> &mask);

/**
 * A convolver object convolves two images
 */
class Convolver {

public:
	virtual ~Convolver();

	/**
	 * Convolves image `src` with the kernel `krn`.
	 * A mask parameter also controls which pixels from the original image
	 * should be convolved. If empty, all pixels are convolved.

	 * @param src The source image
	 * @param src_width The source image's width
	 * @param src_height The source image's height
	 * @param krn The convolution kernel
	 * @param krn_width The kernel's width
	 * @param krn_height The kernel's height
	 * @param mask An optional boolean mask indicating which pixels of the resulting
	 *             image should be convolved
	 * @return The convolved image
	 */
	virtual
	std::vector<double>
	convolve(const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
	         const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
	         const std::vector<bool> &mask) = 0;

};

/**
 * A brute-force convolver.
 */
class BruteForceConvolver : public Convolver {

public:

	std::vector<double> convolve(
	         const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
	         const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
	         const std::vector<bool> &mask) override;

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

	std::vector<double> convolve(
	         const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
	         const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
	         const std::vector<bool> &mask) override;

private:
	std::unique_ptr<FFTPlan> plan;

	std::vector<std::complex<double>> krn_fft;

	bool reuse_krn_fft;
};

#endif /* PROFIT_FFTW */

} /* namespace profit */

#endif /* PROFIT_CONVOLVE_H */