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
#ifndef PROFIT_MOFFAT_H
#define PROFIT_MOFFAT_H

#include "profit/radial.h"

namespace profit
{

/**
 * A Moffat profile
 *
 * The moffat profile has parameters ``fwhm`` and ``con``, and is calculated as
 * follows at radius `r`:
 *
 * @f[
 *   \left[ 1 + \left(\frac{r}{r_d}\right)\right]^{-con}
 * @f]
 *
 * with:
 *
 * @f[
 *   r_d = \frac{fwhm}{2 \sqrt{2^{\frac{1}{con}} - 1}}
 * @f]
 */
class MoffatProfile : public RadialProfile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	MoffatProfile(const Model &model, const std::string &name);

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
	 * Full-width at half maximum of the profiles across the major axis of the
	 * intensity profile.
	 */
	double fwhm;

	/**
	 * %Profile concentration
	 */
	double con;
	// @}

private:
  double fluxfrac(double fraction) const;

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

#endif /* PROFIT_MOFFAT_H */
