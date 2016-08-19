/**
 * Header file for moffat profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham and Rodrigo Tobar
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
#ifndef _MOFFAT_H_
#define _MOFFAT_H_

#include "profit/radial.h"

namespace profit
{

/**
 * A Moffat profile
 *
 * The moffat profile has parameters ``fwhm`` and ``con``, and is calculated as
 * follows for coordinates x/y::
 *
 *   (1+r_factor)^(-c)
 *
 * with::
 *
 *   r_factor = (r/rscale)^2
 *     rscale = fwhm/(2*sqrt( 2^(1/con) - 1))
 *          r = (x^{2+b} + y^{2+b})^{1/(2+b)}
 *          b = box parameter
 */
class MoffatProfile : public RadialProfile {

protected:
	/* All these are inherited from RadialProfile */
	double get_lumtot(double r_box);
	double get_rscale();
	double adjust_acc();
	double adjust_rscale_switch();
	double adjust_rscale_max();
	eval_function_t get_evaluation_function();

public:

	/**
	 * Constructor
	 */
	MoffatProfile();

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/**
	 * Full-width at half maximum of the profiles across the major axis of the
	 * intensity profile.
	 */
	double fwhm;

	/**
	 * Profile concentration
	 */
	double con;

};

} /* namespace profit */

#endif /* _MOFFAT_H_ */
