/**
 * Image convolution implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Rodrigo Tobar
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

#include <stdlib.h>
#include <string.h>

#include "convolve.h"

void profit_normalize(double *image, unsigned int img_width, unsigned int img_height) {

	unsigned int i;
	unsigned int size = img_width * img_height;
	double sum = 0;

	double *in = image;
	for(i=0; i!=size; i++) {
		sum += *in;
		in++;
	}

	in = image;
	for(i=0; i!=size; i++) {
		*in /= sum;
		in++;
	}

}

double *profit_convolve(double *src, unsigned int src_width, unsigned int src_height,
                        double *krn, unsigned int krn_width, unsigned int krn_height,
                        bool *mask, bool replace){

	double pixel;
	unsigned int i, j, k, l;
	unsigned int krn_center_x = krn_width / 2;
	unsigned int krn_center_y = krn_height / 2;
	int src_i, src_j;

	double *convolution = (double *)calloc(src_width * src_height, sizeof(double));

	/* Convolve! */
	/* Loop around the output image first... */
	for (j = 0; j < src_height; j++) {
		for (i = 0; i < src_width; i++) {

			/* Don't convolve this pixel */
			if( mask && !mask[i + j*src_width] ) {
				continue;
			}

			/* ... now loop around the kernel */
			pixel = 0;
			for (l = 0; l < krn_height; l++) {
				for (k = 0; k < krn_width; k++) {

					src_i = (int)i + (int)l - (int)krn_center_x;
					src_j = (int)j + (int)k - (int)krn_center_y;
					if( src_i >= 0 && src_i < src_width &&
					    src_j >= 0 && src_j < src_height ) {
						pixel +=  src[(unsigned int)src_i + (unsigned int)src_j*src_width] * krn[k + l*krn_width];
					}
				}
			}
			convolution[i + j*src_width] = pixel;
		}
	}

	if( replace ) {
		src = memcpy(src, convolution, sizeof(double) * src_width * src_height);
		free(convolution);
		return src;
	}

	return convolution;
}
