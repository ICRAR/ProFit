/**
 * Header file for BrokenExponential profile implementation
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
#ifndef _BROKENEXPONENTIAL_H_
#define _BROKENEXPONENTIAL_H_

#include "profit/radial.h"

namespace profit
{

/**
 * A Broken Exponential profile
 *
 * The Broken Exponential profile has parameters ``h1``, ``h2``, ``rb`` and ``a`` is
 * calculated as follows for coordinates x/y::
 *
 *    inten = exp(-r/h1)*
 *            (1+exp(a*(r-rb)))^((1/a)*(1/h1-1/h2))
 *
 * with::
 *
 *           r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *           B = box parameter
 */
class BrokenExponentialProfile : public RadialProfile {

protected:

	/* All these are inherited from RadialProfile */
	double get_lumtot(double r_box) override;
	double get_rscale() override;
	double adjust_acc() override;
	double adjust_rscale_switch() override;
	double adjust_rscale_max() override;
	eval_function_t get_evaluation_function() override;

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 */
	BrokenExponentialProfile(const Model &model);

	void validate() override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/**
	 * The inner exponential scale length.
	 */
	double h1;

	/**
	 * The outer exponential scale length (must be equal to or less than ``h1``).
	 */
	double h2;

	/**
	 * The break radius.
	 */
	double rb;

	/**
	 * The strength of the truncation as the radius approaches ``rb``.
	 */
	double a;

};

} /* namespace profit */

#endif /* _BROKENEXPONENTIAL_H_ */
