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
typedef double (*eval_function_t)(RadialProfile *, double, double, double, bool);

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
	void initial_calculations();

	/**
	 * Calculates the ``res`` and ``max_rec`` subsampling parameters used for
	 * the first subsampling level of image pixel ``x``/``y``.
	 * Subclasses might want to override this method to provide different
	 * initial subsampling logic.
	 */
	void subsampling_params(double x, double y, unsigned int &res, unsigned int &max_rec);

	/**
	 * Returns the factor by which each resulting image pixel value must be
	 * multiplied to yield the final pixel value. The default implementation
	 * returns the pixel area multiplied by Ie, but subclasses might need to
	 * rescale this.
	 */
	double get_pixel_scale();

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

	RadialProfile();
	void validate();
	void evaluate(double *image);

	/* General parameters */
	double xcen;
	double ycen;
	double mag;
	double ang;
	double axrat;
	double box;

	/*
	 * radius scale, profiles provide it in different ways
	 * via get_rscale()
	 */
	double rscale;

	/* Used to control the subsampling */
	bool rough;
	double acc;
	double rscale_switch;
	unsigned int resolution;
	unsigned int max_recursions;
	bool adjust;

	/* Used to avoid outer regions */
	double rscale_max;

#ifdef PROFIT_DEBUG
	/* record of how many subintegrations we've done */
	std::map<int,int> n_integrations;
#endif

};

} /* namespace profit */

#endif /* _RADIAL_H_ */
