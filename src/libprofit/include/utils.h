/**
 * Header file for utility routines
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

#ifndef _UTILS_H_
#define _UTILS_H_

namespace profit
{

/*
 * Adds the individual values from `src` and `dest` and stores the result
 * in `dest`. Both images must have the same given width and height.
 */
void add_images(double *dest, double *src, unsigned int width, unsigned int height);

/**
 * Normalizes the values of image so their total sum is 1.
 *
 * The values are written back into the image, so if the original needs to be retained
 * then a copy should be supplied.
 */
void normalize(double *image, unsigned int img_width, unsigned int img_height);

} /* namespace profit */

#endif /* _UTILS_H_ */
