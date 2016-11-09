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
 * The CoreSersic profile has parameters `re`, `rb`, `nser`, `a` and `b` and is
 * calculated as follows at radius `r`:
 *
 * @f[
 *    \left[1+\left(\frac{r}{r_b}\right)^{-a}\right]^{\frac{b}{a}}
 *        \exp \left[ -b_n \left( \frac{r^a+{r_b}^a}{{r_e}^a}\right)^{\frac{1}{a n_{ser}}} \right]
 * @f]
 *
 */
class CoreSersicProfile : public RadialProfile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	CoreSersicProfile(const Model &model, const std::string &name);

	void validate() override;

protected:

	/*
	 * ----------------------
	 * Inherited from Profile
	 * ----------------------
	 */
	bool parameter_impl(const std::string &name, double val) override;

	/*
	 * ----------------------------
	 * Inherited from RadialProfile
	 * ----------------------------
	 */
	double get_lumtot(double r_box) override;
	double get_rscale() override;
	double adjust_acc() override;
	double adjust_rscale_switch() override;
	double adjust_rscale_max() override;
	void initial_calculations() override;
	eval_function_t get_evaluation_function() override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/** @name Profile Parameters */
	// @{
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
	// @}

	/* internally calculated when the profile is evaluated */
	double _bn;

private:

	double integrate_at(double r) const;
	double evaluate_at(double x, double y, double r, bool reuse_r) const;

};

} /* namespace profit */

#endif /* _CORESERSIC_H_ */