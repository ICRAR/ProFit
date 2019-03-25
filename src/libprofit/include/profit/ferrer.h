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
#ifndef PROFIT_FERRER_H
#define PROFIT_FERRER_H

#include "profit/config.h"
#include "profit/radial.h"

namespace profit
{

/**
 * A Ferrer profile
 *
 * The ferrer profile has parameters `rout`, `a` and `b` and is
 * calculated as follows at radius `r`:
 *
 * @f[
 *    \left[ 1 - \left(\frac{r}{r_{out}}\right)^{(2-b)} \right]^{a}
 * @f]
 */
class FerrerProfile : public RadialProfile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	FerrerProfile(const Model &model, const std::string &name);

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
	double evaluate_at(double x, double y) const override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/** @name Profile Parameters */
	// @{
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
	// @}

#ifdef PROFIT_OPENCL

protected:
	virtual void add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const override;
	virtual void add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const override;

private:
	template <typename FT>
	void add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const;

#endif /* PROFIT_OPENCL */

};

} /* namespace profit */

#endif /* PROFIT_FERRER_H */
