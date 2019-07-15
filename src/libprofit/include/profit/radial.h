/**
 * Header file for the base radial profile implementation
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
#ifndef PROFIT_RADIAL_H
#define PROFIT_RADIAL_H

#ifdef PROFIT_DEBUG
#include <map>
#endif

#include "profit/config.h"
#include "profit/opencl_impl.h"
#include "profit/profile.h"

namespace profit
{

/**
 * The base class for radial profiles.
 *
 * This class implements the common aspects of all radial profiles, namely:
 *  * High-level evaluation logic
 *  * Region masking
 *  * Translation, rotation, axis ratio and boxing handling
 *  * Pixel subsampling
 *
 * Subclasses are expected to implement a handful of methods that convey
 * profile-specific information, such as the evaluation function for an given
 * x/y profile coordinate and the calculation of the total luminosity of the
 * profile, among others.
 */
class RadialProfile : public Profile {

	friend class FerrerProfile;
	friend class MoffatProfile;
	friend class SersicProfile;

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	RadialProfile(const Model &model, const std::string &name);

	/*
	 * ---------------------------------------------
	 * Pure virtual functions implementations follow
	 * ---------------------------------------------
	 */
	void validate() override;
	void evaluate(Image &image, const Mask &mask, const PixelScale &scale,
	    const Point &offset, double magzero) override;

#ifdef PROFIT_DEBUG
	std::map<int,int> get_integrations();
#endif

protected:

	/**
	 * Calculates a *boxy radius* using the boxiness parameter `box` of this
	 * profile.
	 *
	 * @param x `x` coordinate
	 * @param y `y` coordinate
	 * @return The *boxy* radius `r` for the coordinate (`x`, `y`).
	 */
	double boxy_r(double x, double y) const
	{
		if (box == 0) {
			return std::sqrt(x * x + y * y);
		}
		double box_plus_2 = box + 2.;
		return std::pow(std::pow(std::abs(x), box_plus_2) +
		                    std::pow(std::abs(y), box_plus_2),
		                1. / box_plus_2);
	}

	/**
	 * Calculates the profile value at profile coordinates ``x``/``y``.
	 * @param x The X profile coordinate to evaluate
	 * @param y The Y profile coordinate to evaluate
	 * @return The value of the profile at the given point
	 */
	virtual double evaluate_at(double x, double y) const = 0;

	/**
	 * Performs the initial calculations needed by this profile during the
	 * evaluation phase. Subclasses might want to override this method to add
	 * their own initialization steps.
	 */
	virtual void initial_calculations();

	/**
	 * Calculates the ``res`` and ``max_rec`` subsampling parameters used for
	 * the first subsampling level of image pixel ``x``/``y``.
	 * Subclasses might want to override this method to provide different
	 * initial subsampling logic.
	 */
	virtual void subsampling_params(double x, double y, unsigned int &res, unsigned int &max_rec);

	/**
	 * Returns the factor by which each resulting image pixel value must be
	 * multiplied to yield the final pixel value. The default implementation
	 * returns the pixel area multiplied by Ie, but subclasses might need to
	 * rescale this.
	 *
	 * @param scale the pixel scale as given by the Model when calling @ref evaluate
	 */
	virtual double get_pixel_scale(const PixelScale &scale);

	/*
	 * ------------------------------------------
	 *  Mandatory subclasses methods follow
	 * ------------------------------------------
	 */

	/**
	 * Returns the total luminosity of this profile.
	 */
	virtual double get_lumtot() = 0;

	/**
	 * Returns the value used as ``rscale`` by the radial profile.
	 */
	virtual double get_rscale() = 0;

	/**
	 * Returns an automatically adjusted value for the subsampling accuracy,
	 * which will replace the default or user-given value if users decide to
	 * let the code self-adjust.
	 * The default implementation leaves the accuracy untouched, but subclasses
	 * can override this method.
	 */
	virtual double adjust_acc(double acc);

	/**
	 * Returns an automatically adjusted value for the rscale_switch flag,
	 * which will replace the default or user-given value if users decide to
	 * let the code self-adjust.
	 */
	virtual double adjust_rscale_switch() = 0;

	/**
	 * Returns an automatically adjusted value for the rscale_max flag,
	 * which will replace the default or user-given value if users decide to
	 * let the code self-adjust.
	 */
	virtual double adjust_rscale_max() = 0;

private:
	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/** @name Profile Parameters */
	// @{
	/**
	 * The X center of this profile, in image coordinates
	 */
	double xcen;

	/**
	 * The Y center of this profile, in image coordinates
	 */
	double ycen;

	/**
	 * The magnitude of this profile.
	 */
	double mag;

	/**
	 * The angle by which this profile is rotated. 0 is north, positive is
	 * counterclockwise.
	 */
	double ang;

	/**
	 * The ratio between the two axes, expressed as minor/major.
	 */
	double axrat;

	/**
	 * The *boxiness* of this profile.
	 */
	double box;

	/*
	 * ---------------------------------------
	 * Sub-pixel integration parameters follow
	 * ---------------------------------------
	 */

	/**
	 * Whether perform sub-pixel integration or not.
	 */
	bool rough;

	/**
	 * Target accuracy to achieve during sub-pixel integration
	 */
	double acc;

	/**
	 * Radius (relative to `rscale`) under which sub-pixel integration should
	 * take place
	 */
	double rscale_switch;

	/**
	 * Resolution of the sub-pixel integration: each area to be sub-integrated
	 * is divided in `resolution * resolution` cells.
	 */
	unsigned int resolution;

	/**
	 * Maximum number of recursions that the sub-pixel integration algorithm
	 * should undertake.
	 */
	unsigned int max_recursions;

	/**
	 * Whether this profile should adjust the sub-pixel integration parameters
	 * automatically based on the profile parameters
	 */
	bool adjust;

	/**
	 * Radius (relative to `rscale`) after which the profile is not evaluated
	 * anymore
	 */
	double rscale_max;

	/// Whether the CPU evaluation method should be used, even if an OpenCL
	/// environment has been given (and libprofit has been compiled with OpenCL support)
	bool force_cpu;
	// @}

	/*
	 * radius scale, profiles provide it in different ways
	 * via get_rscale()
	 */
	double rscale;

	/* These are internally calculated at profile evaluation time */
	double _ie;
	double _cos_ang;
	double _sin_ang;
	double _xcen;
	double _ycen;
	double magzero;

	void evaluate_cpu(Image &image, const Mask &mask, const PixelScale &scale);

	void _image_to_profile_coordinates(double x, double y, double &x_prof, double &y_prof);

	double subsample_pixel(double x0, double x1,
	                       double y0, double y1,
	                       unsigned int recur_level,
	                       unsigned int max_recursions,
	                       unsigned int resolution);


#ifdef PROFIT_DEBUG
	/* record of how many subintegrations we've done */
	std::map<int,int> n_integrations;
#endif /* PROFIT_DEBUG */

#ifdef PROFIT_OPENCL

	/**
	 * Indicates whether this profile supports OpenCL evaluation or not
	 * (i.e., implements the required OpenCL kernels)
	 *
	 * @return Whether this profile supports OpenCL evaluation. The default
	 * implementation returns true.
	 */
	virtual bool supports_opencl() const;

	/* Evaluates this radial profile using an OpenCL kernel and floating type FT */
	template <typename FT>
	void evaluate_opencl(Image &image, const Mask &mask, const PixelScale &scale, OpenCLEnvImplPtr &env);

	template <typename FT>
	void add_common_kernel_parameters(unsigned int argIdx, cl::Kernel &kernel) const;

protected:

	/* Add extra parameters to the given kernel, starts with parameter `index` */
	virtual void add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const;
	virtual void add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const;

#endif /* PROFIT_OPENCL */

};

} /* namespace profit */

#endif /* PROFIT_RADIAL_H */