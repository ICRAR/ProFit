/**
 * Header file for King profile implementation
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
#ifndef _KING_H_
#define _KING_H_

#include "profit/radial.h"

namespace profit
{

/**
 * A King profile
 *
 * The King profile has parameters `rc`, `rt` and `a` is
 * calculated as follows at radius `r`:
 *
 * @f[
 *    \left(
 *      \frac{1}{ \left[1 + \left(\frac{r}{r_c}\right)^{2}\right]^{\frac{1}{a}} } -
 *      \frac{1}{ \left[1 + \left(\frac{r_t}{r_c}\right)^{2}\right]^{\frac{1}{a}} }
 *    \right)^{a}
 * @f]
 */
class KingProfile : public RadialProfile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	KingProfile(const Model &model, const std::string &name);

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
	double rc;

	/**
	 * The transition radius of the Sersic profile
	 */
	double rt;

	/**
	 * The power-law of the King.
	 */
	double a;
	// @}

private:

	double integrate_at(double r) const;
	double evaluate_at(double x, double y, double r, bool reuse_r) const;

};

} /* namespace profit */

#endif /* _KING_H_ */