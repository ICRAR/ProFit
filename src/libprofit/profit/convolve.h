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

} /* namespace profit */

#endif /* _CONVOLUTION_H_ */

