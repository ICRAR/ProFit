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
#include <vector>

#include "profit/convolve.h"
#include "profit/exceptions.h"
#include "profit/utils.h"


namespace profit
{

std::vector<double> convolve(
         const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
         const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
         const std::vector<bool> &mask)
{
	return BruteForceConvolver().convolve(src, src_width, src_height, krn, krn_width, krn_height, mask);
}

Convolver::~Convolver()
{
	// no-op
}


std::vector<double> BruteForceConvolver::convolve(
         const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
         const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
         const std::vector<bool> &mask)
{

	double pixel;
	unsigned int i, j, k, l;
	unsigned int krn_half_width = (krn_width - 1) / 2;
	unsigned int krn_half_height = (krn_height - 1) / 2;
	unsigned int krn_size = krn_width * krn_height;
	int src_i, src_j;

	std::vector<double> convolution(src_width * src_height);

	double *out = convolution.data() - 1;
	const double *srcPtr1 = src.data() - 1, *srcPtr2;
	const double *krnPtr;
	std::vector<bool>::const_iterator mask_it = mask.begin();

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
			krnPtr = krn.data() + krn_size - 1;
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

std::vector<double> FFTConvolver::convolve(
        const std::vector<double> &src, unsigned int src_width, unsigned int src_height,
        const std::vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
        const std::vector<bool> &mask)
{

	typedef std::complex<double> complex;

	// Create extended images first
	auto krn_start_x = (src_width - krn_width) / 2;
	auto krn_start_y = (src_height - krn_height) / 2;
	auto ext_width = 2 * src_width;
	auto ext_height = 2 * src_height;
	auto ext_img = extend(src, src_width, src_height, ext_width, ext_height, 0, 0);
	auto ext_krn = extend(krn, krn_width, krn_height, ext_width, ext_height, krn_start_x, krn_start_y);

	// Forward FFTs
	std::vector<complex> src_fft = plan->forward(ext_img);
	if (krn_fft.empty()) {
		krn_fft = plan->forward(ext_krn);
	}

	// element-wise multiplication
	std::transform(src_fft.begin(), src_fft.end(), krn_fft.begin(), src_fft.begin(), std::multiplies<complex>());

	if (!reuse_krn_fft) {
		krn_fft.clear();
	}

	// inverse FFT and scale down
	using namespace std::placeholders;
	std::vector<double> res = plan->backward_real(src_fft);
	std::transform(res.begin(), res.end(), res.begin(),
	               std::bind(std::divides<double>(), _1, res.size()));

	// crop the image to original size
	auto x_offset = src_width / 2;
	auto y_offset = src_height / 2;

	// even image and odd kernel requires slight adjustment
	if (src_width % 2 == 0 and krn_width % 2 == 1) {
		x_offset -= 1;
	}
	if (src_height % 2 == 0 and krn_height % 2 == 1) {
		y_offset -= 1;
	}

	auto cropped = crop(res, ext_width, ext_height, src_width, src_height, x_offset, y_offset);

	// apply mask and good bye
	if (!mask.empty()) {
		std::transform(cropped.begin(), cropped.end(), mask.begin(), cropped.begin(), [](const double i, const bool m) {
			return m ? i : 0.;
		});
	}
	return cropped;
}

#endif /* PROFIT_FFTW */

} /* namespace profit */
