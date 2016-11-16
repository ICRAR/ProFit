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
#ifndef _RADIAL_H_
#define _RADIAL_H_

#ifdef PROFIT_DEBUG
#include <map>
#endif


#include "profit/profit.h"

namespace profit
{

class RadialProfile;

/**
 * The signature that all evaluation functions must follow.
 * Each profile uses a different evaluation function and provides it to the
 * RadialProfile for evaluation.
 *
 * @param profile The profile to evaluate
 * @param x The X profile coordinate to evaluate
 * @param y The Y profile coordinate to evaluate
 * @param r The pre-calculated radius of the profile for coordinate `(x,y)`
 * @param reuse_r Whether the value of `r` is valid (and can be reused) or not.
 * @return The value of the profile at the given point
 */
typedef double (*eval_function_t)(const RadialProfile &profile, double x, double y, double r, bool reuse_r);

/**
 * The case class for radial profiles.
 *
 * This class implements the common aspects of all radial profiles, namely:
 *  * High-level evaluation logic
 *  * Region masking
 *  * Translation, rotation, axis ratio and boxing
 *  * Pixel subsampling
 *
 * Subclasses are expected to implement a handfull of methods that convey
 * profile-specific information, such as the evaluation function for an given
 * x/y profile coordinate and the calculation of the total luminosity of the
 * profile, among others.
 */
class RadialProfile : public Profile {

protected:

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
	 */
	virtual double get_pixel_scale();

	/*
	 * ------------------------------------------
	 *  Mandatory subclasses methods follow
	 * ------------------------------------------
	 */

	/**
	 * Returns the total luminosity of this profile for the given ``r_box``
	 * factor.
	 */
	virtual double get_lumtot(double r_box) = 0;

	/**
	 * Returns the value used as ``rscale`` by the radial profile.
	 */
	virtual double get_rscale() = 0;

	/**
	 * Returns an automatically adjusted value for the subsampling accuracy,
	 * which will replace the default or user-given value if users decide to
	 * let the code self-adjust.
	 */
	virtual double adjust_acc() = 0;

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

	/**
	 * Returns a pointer to the evaluation function, specific to each profile.
	 * The evaluation function takes as input a X/Y coordinate in profile space
	 * and returns a profile value. It receives a pre-computed radius as well
	 * (assuming box == 0) which can be reused to avoid some extra computations.
	 */
	virtual eval_function_t get_evaluation_function() = 0;

	/* These are internally calculated at profile evaluation time */
	double _ie;
	double _cos_ang;
	double _sin_ang;
	eval_function_t _eval_function;

private:

	void _image_to_profile_coordinates(double x, double y, double &x_prof, double &y_prof);

	double subsample_pixel(double x0, double x1,
	                       double y0, double y1,
	                       unsigned int recur_level,
	                       unsigned int max_recursions,
	                       unsigned int resolution);

public:

	/**
	 * Constructor.
	 *
	 * @param model The model this profile belongs to
	 */
	RadialProfile(const Model &);

	/*
	 * ---------------------------------------------
	 * Pure virtual functions implementations follow
	 * ---------------------------------------------
	 */
	void validate() override;
	void evaluate(std::vector<double> &image) override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

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

	/*
	 * radius scale, profiles provide it in different ways
	 * via get_rscale()
	 */
	double rscale;


#ifdef PROFIT_DEBUG
	/* record of how many subintegrations we've done */
	std::map<int,int> n_integrations;
#endif

};

} /* namespace profit */

#endif /* _RADIAL_H_ */
