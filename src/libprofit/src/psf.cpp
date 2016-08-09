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

#include <algorithm>
#include <cmath>

#include "psf.h"
#include "utils.h"

namespace profit
{

void PsfProfile::validate()  {

	if( !this->model->psf ) {
		throw invalid_parameter("No psf present in the model, cannot produce a psf profile");
	}

}

unsigned int bind(double value, unsigned int max) {
	int intval = static_cast<int>(floor(value));
	if( intval < 0 ) {
		return 0;
	}
	unsigned int uintval = static_cast<unsigned int>(intval);
	return std::min(uintval, max);
}

void PsfProfile::evaluate(double *image) {

	int psf_pix_x, psf_pix_y;
	unsigned int i, pix_x, pix_y;
	double x, y, psf_x, psf_y;
	double total = 0;
	double scale = pow(10, -0.4*(this->mag - this->model->magzero));

	/* Making the code more readable */
	double scale_x = model->scale_x;
	double scale_y = model->scale_y;
	double psf_scale_x = model->psf_scale_x;
	double psf_scale_y = model->psf_scale_y;
	unsigned int width = model->width;
	unsigned int height = model->height;
	unsigned int psf_width = model->psf_width;
	unsigned int psf_height = model->psf_height;

	/* Where we start/end applying the psf into the target image */
	double origin_x = this->xcen - psf_width*psf_scale_x/2.;
	double end_x    = this->xcen + psf_width*psf_scale_x/2.;
	double origin_y = this->ycen - psf_height*psf_scale_y/2.;
	double end_y    = this->ycen + psf_height*psf_scale_y/2.;

	/*
	 * We first loop over the pixels of the image, making sure we don't go
	 * outside the image
	 */
	unsigned int x0 = bind(origin_x/scale_x, width - 1);
	unsigned int y0 = bind(origin_y/scale_y, height - 1);
	unsigned int x1 = bind(end_x/scale_x, width - 1);
	unsigned int y1 = bind(end_y/scale_y, height - 1);

	for(pix_y=y0; pix_y <= y1; pix_y++) {

		y = pix_y * scale_y;

		for(pix_x=x0; pix_x <= x1; pix_x++) {

			x = pix_x * scale_x;

			/*
			 * Image pixel (pix_x,pix_y) covers [x:x+scale_x; y:y+scale_y]
			 * We now find out the range of psf pixels that cover the same space
			 */
			int psf_pix_x0 = bind(floor((x - origin_x) / psf_scale_x), psf_width - 1);
			int psf_pix_x1 = bind(floor((x - origin_x + scale_x) / psf_scale_x), psf_height - 1);
			int psf_pix_y0 = bind(floor((y - origin_y) / psf_scale_y), psf_width - 1);
			int psf_pix_y1 = bind(floor((y - origin_y + scale_y) / psf_scale_y), psf_height - 1);

			/* Accumulate the proportional values from the PSF */
			double val = 0;
			for(psf_pix_y = psf_pix_y0; psf_pix_y <= psf_pix_y1; psf_pix_y++) {

				psf_y = psf_pix_y * psf_scale_y + origin_y;

				for(psf_pix_x = psf_pix_x0; psf_pix_x <= psf_pix_x1; psf_pix_x++) {

					psf_x = psf_pix_x * psf_scale_x + origin_x;

					/*
					 * PSF pixel (psf_pix_x,psf_pix_y) covers [psf_x:psf_x+psf_scale_x; psf_y:psf_y+psf_scale_y]
					 *
					 * Now we find the intersection of this PSF pixel area with the
					 * area covered by the image pixel and add its contribution to
					 * the final value of the image pixel.
					 *
					 * On the X coordinate psf_x will always be <= x+scale_x,
					 * and psf_x+psf_scale_x will always be >= x; likewise for the Y
					 * coordinate.
					 *
					 * The contribution will then be given by:
					 */
					double intersect_x = std::min(x + scale_x, psf_x + psf_scale_x) - std::max(x, psf_x);
					double intersect_y = std::min(y + scale_y, psf_y + psf_scale_y) - std::max(y, psf_y);
					val += model->psf[psf_pix_x + psf_pix_y*psf_width] * (intersect_x * intersect_y)/(psf_scale_x * psf_scale_y);

				}
			}

			/* Finally, write down the final value into our pixel */
			image[pix_x + pix_y*width] = val;
			total += val;
		}
	}

	/* We're done applying the ps, now normalize and scale */
	double multiplier = scale;
	if( total != 0 ) {
		multiplier = scale / total;
	}
	double *img_ptr = image;
	for(i=0; i!=height * height; i++, img_ptr++) {
		*img_ptr *= multiplier;
	}

}

PsfProfile::PsfProfile() :
	Profile(),
	xcen(0),
	ycen(0),
	mag(0)
{
	// no-op
}

} /* namespace profit */
