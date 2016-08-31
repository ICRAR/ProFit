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

#include "profit/sky.h"

using namespace std;

namespace profit {

void SkyProfile::validate() {
	/* no-op for the time being, probably check value in range, etc */
	return;
}

void SkyProfile::evaluate(vector<double> &image) {

	/* In case we need to mask some pixels out */
	auto mask_it = model.calcmask.begin();

	/* Fill the image with the background value */
	for(auto &pixel: image) {

		/* Check the calculation mask and avoid pixel if necessary  */
		if( !model.calcmask.empty() ) {
			if( !*mask_it++ ) {
				continue;
			}
		}

		pixel = this->bg;
	}
}

SkyProfile::SkyProfile(const Model &model) :
	Profile(model),
	bg(0.)
{
	// no-op
}

} /* namespace profit */
