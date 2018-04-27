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

#ifndef PROFIT_SKY_H
#define PROFIT_SKY_H

#include "profit/profile.h"

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
	 * @param name The name of this profile
	 */
	SkyProfile(const Model &model, const std::string &name);

	/*
	 * ---------------------------------------------
	 * Pure virtual functions implementations follow
	 * ---------------------------------------------
	 */
	void validate() override;
	void adjust_for_finesampling(unsigned int finesampling) override;
	void evaluate(Image &image, const Mask &mask, const PixelScale &scale, double magzero) override;

protected:

	/*
	 * ----------------------
	 * Inherited from Profile
	 * ----------------------
	 */
	bool parameter_impl(const std::string &name, double value) override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/** @name Profile Parameters */
	// @{
	/**
	 * The value to fill the image with.
	 */
	double bg;
	// @}

	double requested_bg;
};

} /* namespace profit */

#endif /* PROFIT_SKY_H */
