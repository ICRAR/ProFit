/**
 * Header file for CoreSersic profile implementation
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
#ifndef _CORESERSIC_H_
#define _CORESERSIC_H_

#include "profit/radial.h"

namespace profit
{

/**
 * A CoreSersic profile
 *
 * The CoreSersic profile has parameters ``re``, ``rb``, ``nser``, ``a`` and ``b`` and is
 * calculated as follows for coordinates x/y::
 *
 *    (1+(r/rb)^(-a))^(b/a)*
 *        exp(-bn*(((r^a+rb^a)/re^a))^(1/(nser*a)))
 *
 * with::
 *
 *           r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *           B = box parameter
 */
class CoreSersicProfile : public RadialProfile {

protected:

	/* All these are inherited from RadialProfile */
	double get_lumtot(double r_box) override;
	double get_rscale() override;
	double adjust_acc() override;
	double adjust_rscale_switch() override;
	double adjust_rscale_max() override;
	void initial_calculations() override;
	eval_function_t get_evaluation_function() override;

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 */
	CoreSersicProfile(const Model &model);

	void validate() override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/**
	 * The effective radius of the Sersic component
	 */
	double re;

	/**
	 * The transition radius of the Sersic profile
	 */
	double rb;

	/**
	 * The Sersic index of the Sersic profile
	 */
	double nser;

	/**
	 * The strength of transition from inner core to outer Sersic
	 */
	double a;

	/**
	 * The inner power-law of the Core-Sersic.
	 */
	double b;

	/**
	 * The Sersic bn.
	 */
	double _bn;

};

} /* namespace profit */

#endif /* _CORESERSIC_H_ */
