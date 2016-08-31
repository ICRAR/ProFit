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

#include <vector>

namespace profit
{

/**
 * Adds the individual values from `src` and `dest` and stores the result
 * in `dest`. Both images must have the same size.
 */
void add_images(std::vector<double> &dest, const std::vector<double> &src);

/**
 * Normalizes the values of image so their total sum is 1.
 *
 * The values are written back into the image, so if the original needs to be retained
 * then a copy should be supplied.
 */
void normalize(std::vector<double> &image);

/**
 * Computes the quantile of the gamma distribution for ``p`` and ``shape``
 */
double qgamma(double p, double shape);

/**
 * Computes the probability of the gamma distribution for ``q`` and ``shape``
 */
double pgamma(double q, double shape);

/**
 * Computes the gamma function for ``x``
 */
double gammafn(double x);

/**
 * Computes the beta function for ``a`` and ``b``
 */
double beta(double a, double b);

/**
 * A function that can be integrated
 *
 * @param x The domain value used to evaluate the function
 * @param params Additional parameters used to calculate the function
 * @return The value of the integration function at `x`.
 */
typedef double (*integration_func_t)(double x, void *params);

/**
 * Integrates the function `f` on the semi-infinite interval (a, +Inf) using the
 * QAG algorithm (originally from QUADPACK).
 *
 * @param f The function to integrate
 * @param a The beginning of the integration interval
 * @param params A void pointer to any extra data needed by `f`
 * @return The integration result
 */
double integrate_qagi(integration_func_t f, double a, void *params);

/**
 * Integrates the function `f` on the defined interval (a, b) using the
 * QAG algorithm (originally from QUADPACK).
 *
 * @param f The function to integrate
 * @param a The beginning of the integration interval
 * @param b The end of the integration interval
 * @param params A void pointer to any extra data needed by `f`
 * @return The integration result
 */
double integrate_qags(integration_func_t f, double a, double b, void *params);

} /* namespace profit */

#endif /* _UTILS_H_ */
