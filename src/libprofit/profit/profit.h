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

/* M_PI is not part of C/C++, but usually there */
#include <cmath>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/* The override keyword is not supported until gcc 4.7 */
#if defined(__GNUG__) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 7)
#define override
#endif

#include <exception>
#include <string>
#include <vector>

namespace profit
{

/**
 * Exception class thrown when an invalid parameter has been supplied to either
 * a model or a specific profile.
 */
class invalid_parameter : public std::exception
{

public:
	invalid_parameter(const std::string &what);
	invalid_parameter(const invalid_parameter &e);
	~invalid_parameter() throw();
	const char *what() const throw();

private:
	std::string m_what;

};

/* Forward declaration */
class Model;

/**
 * The base profile class
 */
class Profile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 */
	Profile(const Model & model);

	/**
	 * Destructor
	 */
	virtual ~Profile() = 0;

	/**
	 * Performs the initial profile validation, making sure that all parameters
	 * of the profile are correct and can be safely used to create an image.
	 * This function can signal an error by throwing an invalid_parameter exception.
	 */
	virtual void validate() = 0;

	/**
	 * Performs the profile evaluation and saves the resulting image into
	 * the given `image` array. This is the main function of the profile.
	 *
	 * @param image The vector where image values need to be stored.
	 *              Its size is `model.width` * `model.height`.
	 *              The data is organized by rows first, columns later;
	 *               i.e pixel (x,y) is accessed by `image[y*width + x]`
	 */
	virtual void evaluate(std::vector<double> &image) = 0;

	/**
	 * A (constant) reference to the model this profile belongs to
	 */
	const Model &model;

	/**
	 * Whether the resulting image of this profile should be convolved or not.
	 */
	bool convolve;

	/**
	 * The name of this profile
	 */
	std::string name;

};

/**
 * The overall model to be created
 *
 * The model includes the width and height of the image to produce, as well as
 * the resolution to use when performing calculations. Having resolution
 * allows us to specify pixel position with decimal places; e.g., the center
 * point for a given profile.
 */
class Model {

public:


	/**
	 * Constructor
	 *
	 * It creates a new model to which profiles can be added, and that can be
	 * used to calculate an image.
	 */
	Model();

	/**
	 * Destructor.
	 *
	 * It frees all the resources used by the given model, after which it cannot
	 * be used anymore.
	 */
	~Model();

	/**
	 * Creates a new profile for the given name and adds it to the given model.
	 * On success, the new profile is created, added to the model,
	 * and its reference is returned for further customization.
	 * If a profile with the given name is not supported an invalid_parameter
	 * exception is thrown.
	 *
	 * @param profile_name The name of the profile that should be created
	 * @returns A new profile that corresponds to the given name
	 */
	Profile &add_profile(const std::string &profile_name);

	/**
	 * Whether this model contains any profiles or not.
	 *
	 * @return `true` if this module contains at least one profile,
	 * `false` otherwise
	 */
	bool has_profiles() const;

	/**
	 * Calculates an image using the information contained in the model.
	 * The result of the computation is stored in the image field.
	 *
	 * @returns The image created by libprofit. The data is organized by rows
	 *          first, columns later; i.e pixel ``(x,y)`` is accessed by
	 *          ``image[y*width + x]``
	 */
	std::vector<double> evaluate();

	/**
	 * The width of the model to generate
	 */
	unsigned int width;

	/**
	 * The height of the model to generate
	 */
	unsigned int height;

	/**
	 * The X scale; that is, the width of a single pixel in image coordinates
	 */
	double scale_x;

	/**
	 * The Y scale; that is, the height of a single pixel in image coordinates
	 */
	double scale_y;

	/**
	 * The base magnitude applied to all models
	 */
	double magzero;

	/**
	 * The point spread function (psf) to use when convolving images
	 */
	std::vector<double> psf;

	/**
	 * The psf's width
	 */
	unsigned int psf_width;

	/**
	 * The psf's height
	 */
	unsigned int psf_height;

	/**
	 * The PSF's X scale; that is, the width of a single PSF pixel in image
	 * coordinates
	 */
	double psf_scale_x;

	/**
	 * The PSF's Y scale; that is, the height of a single PSF pixel in image
	 * coordinates
	 */
	double psf_scale_y;

	/**
	 * The calculation mask. If given it must be the same size of the expected
	 * output image, and its values are used to limit the profile calculation
	 * only to a given area (i.e., those cells where the value is ``true``).
	 */
	std::vector<bool> calcmask;

private:

	/**
	 * A list of pointers to the individual profiles used to generate the
	 * model's image
	 */
	std::vector<Profile *> profiles;

};

} /* namespace profit */

#endif /* _PROFIT_H_ */
