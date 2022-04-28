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

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/model.h"
#include "profit/psf.h"
#include "profit/utils.h"


namespace profit
{

void PsfProfile::validate()  {

	if( !model.psf ) {
		throw invalid_parameter("No psf present in the model, cannot produce a psf profile");
	}

}

static inline
unsigned int bind(double value, unsigned int max) {
	int intval = static_cast<int>(std::floor(value));
	if( intval < 0 ) {
		return 0;
	}
	auto uintval = static_cast<unsigned int>(intval);
	return std::min(uintval, max);
}

void PsfProfile::evaluate(Image &image, const Mask & /*mask*/, const PixelScale &pixel_scale,
    const Point &offset, double magzero)
{
	using std::floor;
	using std::min;
	using std::max;
	using std::pow;

	double scale = pow(10, -0.4*(this->mag - magzero));

	/* Making the code more readable */
	double scale_x = pixel_scale.x;
	double scale_y = pixel_scale.y;
	double psf_scale_x = model.psf_scale.x;
	double psf_scale_y = model.psf_scale.y;
	unsigned int width = image.getWidth();
	unsigned int height = image.getHeight();
	unsigned int psf_width = model.psf.getWidth();
	unsigned int psf_height = model.psf.getHeight();

	/* Where we start/end applying the psf into the target image */
	double origin_x = this->xcen + offset.x * scale_x - psf_width * psf_scale_x / 2;
	double end_x    = this->xcen + offset.y * scale_y + psf_width * psf_scale_x / 2;
	double origin_y = this->ycen + offset.x * scale_x - psf_height * psf_scale_y / 2;
	double end_y    = this->ycen + offset.y * scale_y + psf_height * psf_scale_y / 2;

	/*
	 * We first loop over the pixels of the image, making sure we don't go
	 * outside the image
	 */
	unsigned int x0 = bind(origin_x/scale_x, width - 1);
	unsigned int y0 = bind(origin_y/scale_y, height - 1);
	unsigned int x1 = bind(end_x/scale_x, width - 1);
	unsigned int y1 = bind(end_y/scale_y, height - 1);

	for(unsigned int pix_y = y0; pix_y <= y1; pix_y++) {

		double y = pix_y * scale_y;

		for(unsigned int pix_x = x0; pix_x <= x1; pix_x++) {

			double x = pix_x * scale_x;

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
			for(int psf_pix_y = psf_pix_y0; psf_pix_y <= psf_pix_y1; psf_pix_y++) {

				double psf_y = psf_pix_y * psf_scale_y + origin_y;

				for(int psf_pix_x = psf_pix_x0; psf_pix_x <= psf_pix_x1; psf_pix_x++) {

					double psf_x = psf_pix_x * psf_scale_x + origin_x;

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
					double intersect_x = min(x + scale_x, psf_x + psf_scale_x) - max(x, psf_x);
					double intersect_y = min(y + scale_y, psf_y + psf_scale_y) - max(y, psf_y);
					val += model.psf[psf_pix_x + psf_pix_y*psf_width] * (intersect_x * intersect_y)/(psf_scale_x * psf_scale_y);

				}
			}

			/* Finally, write down the final value into our pixel */
			image[pix_x + pix_y*width] += val * scale;
		}
	}

}

PsfProfile::PsfProfile(const Model &model, const std::string &name) :
	Profile(model, name),
	xcen(0),
	ycen(0),
	mag(0)
{
	register_parameter("xcen", xcen);
	register_parameter("ycen", ycen);
	register_parameter("mag", mag);
}

} /* namespace profit */
