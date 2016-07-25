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

#include "sky.h"

namespace profit {

void SkyProfile::validate() {
	/* no-op for the time being, probably check value in range, etc */
	return;
}

void SkyProfile::evaluate(double *image) {

	Model *model = this->model;

	/* Setup a pointer to iterate over the calcmask, if any */
	bool *mask_ptr = model->calcmask;
	if( mask_ptr ) {
		mask_ptr -= 1;
	}

	unsigned int i, size = model->width * model->height;

	/* Fill the image with the background value */
	for(i=0; i!=size; i++) {

		/* Check the calculation mask and avoid pixel if necessary  */
		if( model->calcmask ) {
			mask_ptr++;
			if( !*mask_ptr ) {
				continue;
			}
		}

		*image = this->bg;
		image++;
	}
}

SkyProfile::SkyProfile() :
	Profile(),
	bg(0.)
{
	// no-op
}

} /* namespace profit */
