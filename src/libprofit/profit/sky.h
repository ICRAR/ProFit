/**
 * Header fo sky profile implementation
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

#ifndef _SKY_H_
#define _SKY_H_

#include "profit/profit.h"

namespace profit
{

/**
 * A sky profile.
 *
 * This profiles simply fills the image with a constant ``bg`` value, which is
 * given as a parameter.
 */
class SkyProfile : public Profile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 */
	SkyProfile(const Model &model);

	/*
	 * ---------------------------------------------
	 * Pure virtual functions implementations follow
	 * ---------------------------------------------
	 */
	void validate() override;
	void evaluate(std::vector<double> &image) override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/**
	 * The value to fill the image with.
	 */
	double bg;

};

} /* namespace profit */

#endif /* _SKY_H_ */
