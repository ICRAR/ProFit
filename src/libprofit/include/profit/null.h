/**
 * Header file for null profile implementation
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
#ifndef PROFIT_NULL_H
#define PROFIT_NULL_H

#include "profit/profile.h"

namespace profit
{

/**
 * A null profile.
 *
 * The null profiles has no parameters, and leaves the incoming input
 * image untouched. It is only useful for testing purposes.
 */
class NullProfile : public Profile {
public:
	NullProfile(const Model &model, const std::string &name) : Profile(model, name) {}
	void validate() override {};
	void evaluate(Image &image, const Mask &mask, const PixelScale &scale, double magzero) {};
};

} /* namespace profit */

#endif /* PROFIT_NULL_H */
