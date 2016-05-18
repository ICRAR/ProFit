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

#include <stdlib.h>

#include "sky.h"

static
void profit_init_sky(profit_profile *profile, profit_model *model) {
	/* no-op for the time being, probably check value in range, etc */
	return;
}

static
void profit_make_sky(profit_profile *profile, profit_model *model, double *image) {
	unsigned int i, j;
	profit_sky_profile *sky_p = (profit_sky_profile *)profile;

	/* Fill the image with the background value */
	for(i=0; i!=model->width; i++) {
		for(j=0; j!=model->height; j++) {
			image[j*model->width + i] = sky_p->bg;
		}
	}
}

profit_profile *profit_create_sky() {
	profit_sky_profile *p = (profit_sky_profile *)malloc(sizeof(profit_sky_profile));
	p->profile.init_profile = &profit_init_sky;
	p->profile.make_profile = &profit_make_sky;
	p->bg = 0.;
	return (profit_profile *)p;
}
