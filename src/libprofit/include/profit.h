/**
 * Header file with main libprofit structures and functions
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

#ifndef _PROFIT_H_
#define _PROFIT_H_

#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

struct _profit_model;

/**
 * The base profile structure
 *
 * Specific profile structures *must* declare a profit_profile structure as its
 * first element, so the resulting memory can be addressed both as a generic
 * profit_profile or as the specific profile structure.
 */
typedef struct _profit_profile {

	/**
	 * The name of this profile
	 */
	const char *name;

	/**
	 * Whether the resulting image of this profile should be convolved or not.
	 */
	bool convolve;

	/**
	 * A pointer to the function that performs the initial profile validation,
	 * making sure that all parameters of the profile are correct and can be
	 * safely used to create an image.
	 * This function can signal an error by setting a value in the error member
	 * of this structure.
	 */
	void (* validate_profile)(struct _profit_profile *profile, struct _profit_model *model);

	/**
	 * A pointer to the function that performs the profile evaluation.
	 * This is the main function of the profile.
	 */
	void (* evaluate_profile)(struct _profit_profile *profile, struct _profit_model *model, double *image);

	/**
	 * An error string indicating that an error related to this profile was
	 * detected. The error string can be set either during the profile
	 * initialization or during the image creation process. Users should check
	 * that there is no error in any of the profiles after making a model.
	 */
	char *error;

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
	 * The point spread function (psf) to use when convolving images
	 */
	double *psf;

	/**
	 * The psf's width
	 */
	unsigned int psf_width;

	/**
	 * The psf's height
	 */
	unsigned int psf_height;

	/*
	 * Used to limit the profile calculation only to a given area
	 */
	bool *calcmask;

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

	/**
	 * An error string indicating that there is something wrong with the model.
	 * Users should check that there is no error after making a model.
	 */
	char *error;

} profit_model;

/**
 * Creates a new model to which profiles can be added, and that can
 * be used to calculate an image.
 */
profit_model *profit_create_model(void);

/**
 * Creates a new profile for the given name and adds it to the given model.
 * On success, the new profile is created, added to the model,
 * and its reference is returned for further customization.
 * On failure (i.e., if a profile with the given name is not supported) NULL is
 * returned and no profile is added to the model.
 */
profit_profile *profit_create_profile(profit_model *model, const char *profile_name);

/**
 * Calculates an image using the information contained in the model.
 * The result of the computation is stored in the image field.
 */
void profit_eval_model(profit_model *model);

/**
 * Returns the first error string found either on the model itself or in any of
 * it profiles. This method should be called on the model right after invoking
 * profit_eval_model to make sure that no errors were found during the process.
 * If NULL is returned it means that no errors were found and that the image
 * stored in the model is valid.
 */
char *profit_get_error(profit_model *model);

/**
 * Frees all the resources used by the given model, after which it cannot be
 * used anymore.
 */
void profit_cleanup(profit_model *model);

#ifdef __cplusplus
}
#endif

#endif /* _PROFIT_H_ */
