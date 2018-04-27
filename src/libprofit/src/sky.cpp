/**
 * Sky profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Rodrigo Tobar
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

#include <vector>

#include "profit/common.h"
#include "profit/model.h"
#include "profit/sky.h"


namespace profit {

void SkyProfile::validate() {
	/* no-op for the time being, probably check value in range, etc */
	return;
}


void SkyProfile::adjust_for_finesampling(unsigned int finesampling)
{
	bg = requested_bg / (finesampling * finesampling);
}

void SkyProfile::evaluate(Image &image, const Mask &mask, const PixelScale &scale, double magzero) {

	/* In case we need to mask some pixels out */
	auto mask_it = mask.begin();

	/* Fill the image with the background value */
	for(auto &pixel: image) {

		/* Check the calculation mask and avoid pixel if necessary  */
		if( mask && !*mask_it++ ) {
			continue;
		}

		pixel = this->bg;
	}
}

SkyProfile::SkyProfile(const Model &model, const std::string &name) :
	Profile(model, name),
	bg(0.),
	requested_bg(0.)
{
	// no-op
}

bool SkyProfile::parameter_impl(const std::string &name, double val) {

	if( Profile::parameter_impl(name, val) ) {
		return true;
	}

	if( name == "bg" ) {
		this->requested_bg = val;
		return true;
	}

	return false;
}

} /* namespace profit */
