/**
 * OpenMP utilities for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2018
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

#ifndef PROFIT_OMP_UTILS_H_
#define PROFIT_OMP_UTILS_H_

#include "profit/common.h"

namespace profit {

/**
 * Runs @p f over each point ``(i, j)`` of a grid of width @p width and height
 * @p height using @p threads OpenMP threads. If no OpenMP support is found, @p f
 * is called sequentially in row-first order (that is, ``i`` values changing
 * more rapidly than ``j`` values).
 *
 * @param threads The number of OpenMP threads to use
 * @param width The width of the grid of points
 * @param height The height of the grid of points
 * @param f The function to evaluate on each point. It should receive ``i`` and
 * ``j`` as arguments
 */
template <typename Callable>
void omp_2d_for(int threads, unsigned int width, unsigned int height, Callable &&f)
{
	// Don't even try using OpenMP if serial execution was requested
	// Setting up OpenMP incurs into some costs that can be avoided with this
	if (threads <= 1) {
		for (unsigned int j = 0; j < height; j++) {
			for (unsigned int i = 0; i < width; i++) {
				f(i, j);
			}
		}
		return;
	}
#if _OPENMP >= 200805 // OpenMP 3.0
#pragma omp parallel for collapse(2) schedule(dynamic, 10) if(threads > 1) num_threads(threads)
	for (unsigned int j = 0; j < height; j++) {
		for (unsigned int i = 0; i < width; i++) {
			f(i, j);
		}
	}
#elif _OPENMP >= 200203 // OpenMP 2.0. No "collapse", signed int loop variable
#pragma omp parallel for schedule(dynamic, 10) if(threads > 1) num_threads(threads)
	for (int x = 0; x < int(width * height); x++) {
		unsigned int i = x % width;
		unsigned int j = x / width;
		f(i, j);
	}
#else
	UNUSED(threads);
	for (unsigned int j = 0; j < height; j++) {
		for (unsigned int i = 0; i < width; i++) {
			f(i, j);
		}
	}
#endif // _OPENMP
}

}  // namespace profit

#endif /* PROFIT_OMP_UTILS_H_ */
