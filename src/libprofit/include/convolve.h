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
#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_

#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Convolves image src with the kernel krn.
 *
 * Both the source image and the kernel need to specify their width and height.
 * Depending on the value of the replace parameter, the same src image will be
 * used to store the convolution result (if replace != 0), or a new vector will
 * be allocated and filled instead.
 */
double *profit_convolve(double *src, unsigned int src_width, unsigned int src_height,
                        double *krn, unsigned int krn_width, unsigned int krn_height,
                        bool replace);

#ifdef __cplusplus
}
#endif

#endif /* _CONVOLUTION_H_ */

