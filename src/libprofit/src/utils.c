/**
 * Utility routines for libprofit
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

void profit_add_images(double *dest, double *src,
                       unsigned int width, unsigned int height) {

	for(unsigned int i=0; i != width*height; i++, dest++, src++) {
		*dest += *src;
	}

}

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
