/**
 * Header file for ferrer profile implementation
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
#ifndef _FERRER_H_
#define _FERRER_H_

#include "profit/radial.h"

namespace profit
{

/**
 * A Ferrer profile
 *
 * The ferrer profile has parameters ``rout``, ``a`` and ``b`` and is
 * calculated as follows for coordinates x/y::
 *
 *    (1-r_factor)^(a)
 *
 * with::
 *
 *    r_factor = (r/rout)^(2-b)
 *           r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *           B = box parameter
 */
class FerrerProfile : public RadialProfile {

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
	FerrerProfile();

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/**
	 * The outer truncation radius
	 */
	double rout;

	/**
	 * The global power-law slope to the profile center
	 */
	double a;

	/**
	 * The strength of the truncation as the radius approaches ``rout``.
	 */
	double b;

};

} /* namespace profit */

#endif /* _FERRER_H_ */
