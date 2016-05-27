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
#include <string.h>
#include <stdlib.h>

#include "psf.h"

static
void profit_init_psf(profit_profile *profile, profit_model *model)  {
	profit_psf_profile *psf = (profit_psf_profile *)profile;

	if( !model->psf ) {
		psf->profile.error = strdup("No psf present in the model, cannot produce a psf profile");
		return;
	}
	psf->scale = pow(10, -0.4*(psf->mag - model->magzero));

}

static inline
void profit_psf_normalize_and_apply(profit_psf_profile *psf, profit_model *model, double *image,
                                    double *psf_img, unsigned int psf_w, unsigned int psf_h,
												int target_x, int target_y) {

	unsigned int i, j, img_x, img_y;

	double total = 0;
	for(j=0; j!=psf_h; j++) {
		for(i=0; i!=psf_w; i++) {
			total += psf_img[i + j*psf_w];
		}
	}

	double scale = psf->scale / total;
	for(j=0; j!=psf_h; j++) {

		/* Don't draw outside the boundaries of the full image */
		if( (int)j+target_y < 0 ) {
			continue;
		}
		img_y = j + (unsigned int)target_y;
		if( img_y >= model->height ) {
			break;
		}

		for(i=0; i!=psf_w; i++) {

			/* Don't draw outside the boundaries of the full image */
			if( (int)i+target_x < 0 ) {
				continue;
			}
			img_x = i + (unsigned int)target_x;
			if( img_x >= model->width ) {
				break;
			}

			image[img_x + img_y*model->width] = psf_img[i + j*psf_w] * scale;
		}
	}

}

static
void profit_make_psf(profit_profile *profile, profit_model *model, double *image) {

	/*
	 * TODO: This method still doesn't take into account the image xbin/ybin
	 */

	unsigned int i, j;
	profit_psf_profile *psf = (profit_psf_profile *)profile;

	/*
	 * The PSF is not simply put "as is" in the nearest position of the desired
	 * center. Its values are interpolated instead to take into account the
	 * sub-pixel distance needed to reach the center.
	 *
	 * We avoid this extra calculation only if the target position of the psf's
	 * origin is located exactly on a pixel crossing, and if the psf's
	 * dimensions are both even. This would ensure that each pixel of the PSF
	 * corresponds exactly to one pixel on the target image, allowing us to have
	 * a direct copy of values.
	 */
	double psf_origin_x = psf->xcen - model->psf_width/2.;
	double psf_origin_y = psf->ycen - model->psf_height/2.;
	if( (model->psf_width % 2 == 0 && model->psf_height % 2 == 0) && \
       (floor(psf_origin_x) == psf_origin_x || ceil(psf_origin_x) == psf_origin_x) && \
	    (floor(psf_origin_y) == psf_origin_y || ceil(psf_origin_y) == psf_origin_y) ) {

		profit_psf_normalize_and_apply(psf, model, image,
		                               model->psf, model->psf_width, model->psf_height,
		                               (int)psf_origin_x, (int)psf_origin_y);

		return;
	}

	/*
	 * Each image pixel is now divided into four regions because of its
	 * intersection with the PSF pixels:
	 *
	 *             |
	 *    ---------|------
	 *    |        |      |
	 *    |   a3   |  a4  |  yd2
	 *  ___________|________
	 *    |        |      |
	 *    |        |      |
	 *    |   a1   |  a2  |  yd1
	 *    |        |      |
	 *    ---------|-------
	 *        xd1  |  xd2
	 *
	 * We average the four areas to obtain the value of the pixel. The areas are
	 * all the same on each image pixel so we calculate them once.
	 */

	double xd1 = psf->xcen - floor(psf->xcen);
	double xd2 = 1 - xd1;
	double yd1 = psf->ycen - floor(psf->ycen);
	double yd2 = 1 - yd1;
	double a1 = xd1 * yd1;
	double a2 = xd2 * yd2;
	double a3 = xd1 * yd2;
	double a4 = xd2 * yd2;

	unsigned int new_psf_w = model->psf_width + 1;
	unsigned int new_psf_h = model->psf_height + 1;
	double *new_psf = (double *)calloc(new_psf_w * new_psf_h, sizeof(double));

	for(j=0; j!=new_psf_h; j++) {
		for(i=0; i!=new_psf_w; i++) {

			/* The borders of the target image area use less psf pixels */
			double psf_val = 0;
			if( i != 0 && j != 0 ) {
				psf_val += model->psf[i-1 + (j-1)*model->psf_width] * a1;
			}
			if( i != 0 && j != (new_psf_h - 1) ) {
				psf_val += model->psf[i-1 + j*model->psf_width] * a3;
			}
			if( i != (new_psf_w - 1) && j != 0 ) {
				psf_val += model->psf[i + (j-1)*model->psf_width] * a2;
			}
			if( i != (new_psf_w - 1) && j != (new_psf_h - 1) ) {
				psf_val += model->psf[i + j*model->psf_width] * a4;
			}

			new_psf[i + j*new_psf_w] = psf_val;
		}
	}

	profit_psf_normalize_and_apply(psf, model, image,
                                  new_psf, new_psf_w, new_psf_h,
	                               (int)floor(psf_origin_x), (int)floor(psf_origin_y));

	free(new_psf);

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
