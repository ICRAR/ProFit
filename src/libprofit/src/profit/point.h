/**
 * Point classes definitions
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2021
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


#ifndef PROFIT_POINT_H
#define PROFIT_POINT_H

#include "profit/coordinates.h"

namespace profit {

/// A point in a discrete 2-dimensional surface
class Point : public discrete_2d_coordinate
{
public:
	using discrete_2d_coordinate::discrete_2d_coordinate;
};

/// A point in a continuous 2-dimensional surface
class ContPoint : public continuous_2d_coordinate<ContPoint>
{
public:
	using continuous_2d_coordinate<ContPoint>::continuous_2d_coordinate;
};

}  // namespace profit

#endif /* PROFIT_POINT_H */