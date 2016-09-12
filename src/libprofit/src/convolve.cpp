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

#include <vector>

using namespace std;

namespace profit
{

vector<double>
convolve(const vector<double> &src, unsigned int src_width, unsigned int src_height,
         const vector<double> &krn, unsigned int krn_width, unsigned int krn_height,
         const vector<bool> &mask){

	double pixel;
	unsigned int i, j, k, l;
	unsigned int krn_half_width = (krn_width - 1) / 2;
	unsigned int krn_half_height = (krn_height - 1) / 2;
	unsigned int krn_size = krn_width * krn_height;
	int src_i, src_j;

	vector<double> convolution(src_width * src_height);

	double *out = convolution.data() - 1;
	const double *srcPtr1 = src.data() - 1, *srcPtr2;
	const double *krnPtr;
	vector<bool>::const_iterator mask_it = mask.begin();

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

} /* namespace profit */
