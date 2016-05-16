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

#include "convolve.h"

double *profit_convolve(double *src, unsigned int src_width, unsigned int src_height,
                        double *krn, unsigned int krn_width, unsigned int krn_height,
                        bool replace){

	double pixel;
	unsigned int i, j, k, l;
	unsigned int extrax = krn_width % 2;
	unsigned int extray = krn_height % 2;

	unsigned int conv_width = src_width + krn_width - extrax;
	unsigned int conv_height = src_height + krn_height - extray;
	double *convolution = (double *)calloc(conv_width * conv_height, sizeof(double));

	/* Convolve! */
	/* Loop around the image first... */
	for (j = 0; j < src_height; j++) {
		for (i = 0; i < src_width; i++) {
			pixel = src[i + j*src_width];

			/* ... now loop around the kernel */
			for (l = 0; l < krn_height; l++) {
				for (k = 0; k < krn_width; k++) {
					convolution[k+i + (l + j)*conv_width] += pixel * krn[k + l*krn_width];
				}
			}
		}
	}

	/* Cut out to the original size and return */
	double *cutout = src;
	if( !replace ) {
		cutout = (double *)malloc(sizeof(double) * src_width * src_height);
	}
	for (j = 0; j < src_height; j++) {
		for (i = 0; i < src_width; i++) {
			src[i + j*src_width] = convolution[i+(krn_width-extrax)/2 + (j+(krn_height-extray)/2) * conv_width];
		}
	}
	free(convolution);

	return cutout;
}
