/**
 * Header file with common definitions for libprofit
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

#ifndef PROFIT_COMMON_H
#define PROFIT_COMMON_H

/* M_PI is not part of C/C++, but usually there */
#include <cmath>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

/* The override keyword is not supported until gcc 4.7 */
#if defined(__GNUG__) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 7)
# define override
#endif

namespace profit {
	typedef unsigned int nsecs_t;
}

#endif /* PROFIT_COMMON_H */
