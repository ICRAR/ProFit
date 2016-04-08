/**
 * Header file with main libprofit structures and functions
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2014
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

#ifndef _PROFIT_H_
#define _PROFIT_H_

#ifdef __cplusplus
extern "C"
{
#endif

struct _profit_model;

typedef struct _profit_profile {
	char *name;
	int (* init_profile)(struct _profit_profile *profile, struct _profit_model *model);
	void (* make_profile)(struct _profit_profile *profile, struct _profit_model *model, double *image);
} profit_profile;

/**
 * The overall model to be created
 *
 * The model includes the width and height of the image to produce, as well as
 * the resolution to use when performing calculations. Having resolution
 * allows us to specify pixel position with decimal places; e.g., the center
 * point for a given profile.
 */
typedef struct _profit_model {

	/**
	 * The width of the model to generate
	 */
	unsigned int width;

	/**
	 * The height of the model to generate
	 */
	unsigned int height;

	/**
	 * The horizontal resolution to use when generating the model
	 */
	unsigned int res_x;

	/**
	 * The vertical resolution to use when generating the model
	 */
	unsigned int res_y;

	/* These are calculated from the widht/height and res fields */
	double xbin;
	double ybin;

	double magzero;

	/**
	 * The image created by libprofit.
	 *
	 * The data is organized by rows first, columns later;
	 * i.e pixel (x,y) is accessed by image[y*width + x]
	 */
	double *image;

	/**
	 * The number of profiles used to generate the model's image
	 */
	unsigned int n_profiles;

	/**
	 * A list of pointers to the individual profiles used to generate the
	 * model's image
	 */
	profit_profile **profiles;

} profit_model;

/**
 * Main entry point routine. It calculates an image using the parameters
 * contained in model (see the definition of the structure for more details).
 * The result of the computation is stored in the image field.
 */
int profit_make_model(profit_model *model);

/**
 * Gets a new profile structure for the given name. If a profile with the given
 * name is not supported then NULL is returned.
 */
profit_profile *profit_get_profile(const char *name);

/**
 * Returns a new model with all n profiles given as arguments bound to it.
 */
profit_model *profit_get_model(unsigned int n, ...);

#ifdef __cplusplus
}
#endif

#endif /* _PROFIT_H_ */
