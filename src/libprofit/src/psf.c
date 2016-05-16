/**
 * PSF profile implementation
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

#include <math.h>
#include <stdlib.h>

#include "psf.h"

static
void profit_init_psf(profit_profile *profile, profit_model *model)  {
	profit_psf_profile *psf = (profit_psf_profile *)profile;
	psf->scale = pow(10, -0.4*(psf->mag - model->magzero));
}

static
void profit_make_psf(profit_profile *profile, profit_model *model, double *image) {

	unsigned int i, j, image_x, image_y;
	profit_psf_profile *psf = (profit_psf_profile *)profile;

	unsigned int xcen_pixel = (unsigned int)round(psf->xcen/model->xbin);
	unsigned int ycen_pixel = (unsigned int)round(psf->ycen/model->ybin);
	unsigned int i_0 = xcen_pixel - model->psf_width/2;
	unsigned int j_0 = ycen_pixel - model->psf_height/2;

	for(i=0; i < model->psf_width; i++) {
		image_x = i + i_0;
		if( image_x > model->width ) {
			break;
		}
		for(j=0; j < model->psf_height ; j++) {
			image_y = j + j_0;
			if( image_y > model->height ) {
				break;
			}
			image[image_y*model->width + image_x] = model->psf[i + j*model->psf_width] * psf->scale;
		}
	}
}

profit_profile *profit_create_psf() {

	profit_psf_profile *psf = (profit_psf_profile *)malloc(sizeof(profit_psf_profile));
	psf->profile.init_profile = &profit_init_psf;
	psf->profile.make_profile = &profit_make_psf;

	psf->xcen = 0;
	psf->ycen = 0;
	psf->mag  = 0;
	return (profit_profile *)psf;
}
