/**
 * Header file for the main Model class of libprofit
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

#ifndef PROFIT_MODEL_H
#define PROFIT_MODEL_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "profit/config.h"
#include "profit/common.h"
#include "profit/convolve.h"
#include "profit/opencl.h"
#include "profit/profile.h"

namespace profit
{

/* Forward declaration */
struct ProfileStats;

/**
 * A pair of double values that specify a the scale of a pixel in a particular image
 * with respect to that image's coordinate space.
 */
typedef std::pair<double, double> PixelScale;

/**
 * The overall model to be created
 *
 * The model includes the width and height of the image to produce, as well as
 * the resolution to use when performing calculations. Having resolution
 * allows us to specify pixel position with decimal places; e.g., the center
 * point for a given profile.
 */
class PROFIT_API Model {

public:

	/**
	 * Constructor
	 *
	 * It creates a new model to which profiles can be added, and that can be
	 * used to calculate an image.
	 */
	Model(unsigned int width = 0, unsigned int height = 0);

	/**
	 * Creates a new profile for the given name and adds it to the given model.
	 * On success, the new profile is created, added to the model,
	 * and its reference is returned for further customization.
	 * If a profile with the given name is not supported an invalid_parameter
	 * exception is thrown.
	 *
	 * @param profile_name The name of the profile that should be created
	 * @returns A pointer to the new profile that corresponds to the given name
	 */
	ProfilePtr add_profile(const std::string &profile_name);

	/**
	 * Whether this model contains any profiles or not.
	 *
	 * @return `true` if this module contains at least one profile,
	 * `false` otherwise
	 */
	bool has_profiles() const;

	/**
	 * Calculates an image using the information contained in the model.
	 * The result of the computation is returned as an Image, which may be of a
	 * different size from the one originally requested if the user set this
	 * model's ``crop`` property to ``false`` (via @ref set_crop). If users want to know the offset
	 * at which the image resulting of evaluating this Model with its configured
	 * parameters is with respect to the Image value returned by this method,
	 * hen they must provide a Point in @p offset_out, which will contain
	 * the information after the method returns.
	 *
	 * In other words, the Image returned by this method can be bigger than the
	 * Model's dimensions if the user requested this Model to return a non-cropped Image.
	 *
	 * @param offset_out The potential offset with respect to the image returned by
	 * this method at which the image of this Model's dimensions can be found.
	 *
	 * @returns The image created by libprofit.
	 */
	Image evaluate(Point &offset_out = NO_OFFSET);

#ifdef PROFIT_DEBUG
	std::map<std::string, std::map<int, int>> get_profile_integrations() const;
#endif

	/**
	 * Return a map of all profile statistics.
	 *
	 * @return A map indexed by profile name with runtime statistics
	 */
	std::map<std::string, std::shared_ptr<ProfileStats>> get_stats() const;

	/**
	 * Sets the dimensions of the model image to generate
	 * @param dimensions The dimensions of the model image to generate
	 */
	void set_dimensions(const Dimensions &dimensions) {
		requested_dimensions = dimensions;
	}

	/**
	 * Sets the finesampling factor to use in this Model
	 */
	void set_finesampling(unsigned int finesampling) {
		if (finesampling == 0) {
			throw std::invalid_argument("finesampling must be > 0");
		}
		this->finesampling = finesampling;
	}

	/**
	 * Sets the PSF image that this Model should use
	 * @param psf The PSF image that this Model should use
	 */
	void set_psf(const Image &psf) {
		this->psf = psf;
	}

	/**
	 * @see set_psf(const Image &psf)
	 */
	void set_psf(Image &&psf) {
		this->psf = psf;
	}

	/**
	 * Sets the pixel scale of the generated model image.
	 *
	 * The image scale is the width (and height) of a single
	 * pixel in image coordinates.
	 *
	 * @param scale The pixel scale of the model image
	 */
	void set_image_pixel_scale(const PixelScale &scale) {
		this->scale = scale;
	}

	/**
	 * Returns the pixel scale of the generated model image.
	 *
	 * @return the image scale of the generated model image
	 * @see set_image_pixel_scale(double, double)
	 */
	PixelScale get_image_pixel_scale() const {
		return scale;
	}

	/**
	 * Sets the PSF's pixel scale.
	 *
	 * @param scale The pixel scale of the PSF
	 * @see set_image_pixel_scale(double, double)
	 */
	void set_psf_pixel_scale(const PixelScale &scale) {
		psf_scale = scale;
	}

	/**
	 * Returns the pixel scale of the PSF.
	 *
	 * @return the image scale of the generated model image
	 * @see set_psf_pixel_scale(double, double)
	 */
	PixelScale get_psf_pixel_scale() const {
		return psf_scale;
	}


	/**
	 * Sets the base magnitude to be applied to all profiles.
	 *
	 * @param magzero The base magnitude to be applied to all profiles.
	 */
	void set_magzero(double magzero) {
		this->magzero = magzero;
	}

	/**
	 * Set the calculation mask. If given it must be the same size of the expected
	 * output image, and its values are used to limit the profile calculation
	 * only to a given area (i.e., those cells where the value is ``true``).
	 *
	 * @param mask The mask to use to limit profile calculations
	 */
	void set_mask(const Mask &mask) {
		this->mask = mask;
	}

	/**
	 * @see set_mask(const Mask &mask)
	 */
	void set_mask(Mask &&mask) {
		this->mask = mask;
	}

	/**
	 * Set a convolver for this Model.
	 * A convolver is an object used to carry out the convolution, if necessary.
	 * If a convolver is present before calling `evaluate` then it is used.
	 * If missing and one is required, a new one is created internally.
	 */
	void set_convolver(ConvolverPtr convolver) {
		this->convolver = convolver;
	}

	/**
	 * Set the cropping flag.
	 *
	 * Due to their internal workings, some convolvers produce actually bigger
	 * which are (by default) cropped to the size of the original images created
	 * by the profiles. If this option is set to true, then the result of the
	 * convolution will *not* be cropped, meaning that the result of the model
	 * evaluation will be bigger than what was originally requested.
	 *
	 * @param crop Whether this model returns a cropped image (default) or not
	 * to the user.
	 */
	void set_crop(bool crop) {
		this->crop = crop;
	}

	/**
	 * Sets the dry run flag. The dry run flag determines whether the actual
	 * evaluation of profiles should be skipped or not; if skipped profile validation
	 * still takes place.
	 *
	 * @param dry_run Whether evaluation of profiles should take place (default)
	 * or not
	 */
	void set_dry_run(bool dry_run) {
		this->dry_run = dry_run;
	}

	/**
	 * Set the return finesampled flag.
	 *
	 * When users set a finesampling factor on the model (via
	 * set_finesampling()) the image calculated by this model will have bigger
	 * dimensions than those originally set in the Model. This flag controls
	 * whether this bigger image should be returned (default behavior), or
	 * whether a smaller version of the image with dimensions equals to the ones
	 * requested (plus any padding introduced by convolution) should be
	 * returned. If a smaller image is returned, the total flux of the image is
	 * still preserved.
	 *
	 * @param return_finesampled Whether this model should return finesampled
	 * images as-is (`true`) or if they should be downsampled to match the
	 * original model dimensions.
	 */
	void set_return_finesampled(bool return_finesampled) {
		this->return_finesampled = return_finesampled;
	}

	void set_opencl_env(const OpenCLEnvPtr &opencl_env) {
		this->opencl_env = opencl_env;
	}

	OpenCLEnvPtr get_opencl_env() const {
		return opencl_env;
	}

	/**
	 * Sets the maximum number of OpenMP threads to use to evaluate the profiles
	 * contained in this model. 0 threads means that no OpenMP support
	 * has been requested.
	 *
	 * @param omp_threads the number of OpenMP threads to use for profile evaluation
	 */
	void set_omp_threads(unsigned int omp_threads) {
		this->omp_threads = omp_threads;
	}

	/**
	 * Returns the number of OpenMP threads this Model has been configured to work with
	 * @return the number of OpenMP threads this Model has been configured to work with
	 */
	unsigned int get_omp_threads() {
		return this->omp_threads;
	}

	/**
	 * The Point object that indicates that users don't want to retrieve back
	 * the potential image offset when calling @ref evaluate(Point &)
	 */
	static Point NO_OFFSET;

private:

	Dimensions requested_dimensions;
	unsigned int finesampling;
	PixelScale scale;
	double magzero;
	Image psf;
	PixelScale psf_scale;
	Mask mask;
	ConvolverPtr convolver;
	bool crop;
	bool dry_run;
	bool return_finesampled;
	OpenCLEnvPtr opencl_env;
	unsigned int omp_threads;
	std::vector<ProfilePtr> profiles;

	friend class PsfProfile;
	friend class RadialProfile;
};

} /* namespace profit */

#endif /* PROFIT_MODEL_H */
