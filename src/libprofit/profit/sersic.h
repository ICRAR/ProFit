/**
 * Header file for sersic profile implementation
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
#ifndef PROFIT_SERSIC_H
#define PROFIT_SERSIC_H

#include "profit/radial.h"

namespace profit
{

/**
 * A Sersic profile
 *
 * The sersic profile has parameters `nser` and `re` and is calculated as
 * follows at radius `r`:
 *
 * @f[
 *   \exp \left\{ -b_n \left[ \left(\frac{r}{r_e}\right)^{\frac{1}{n_{ser}}} - 1 \right] \right\}
 * @f]
 */
class SersicProfile : public RadialProfile {

public:

	/*
	 * An enumeration of (inverse) power exponents used internally when
	 * evaluating the sersic profile.
	 * This enumeration allows us to write optimized versions of the evaluation
	 * function particular nser/box combinations, avoiding to call pow() and
	 * calling other functions instead (e.g., sqrt which in x86_64 is
	 * implemented in hardware).
	 */
	enum rfactor_invexp_t {
		general,
		pointfive,
		one,
		two,
		three,
		four,
		eight,
		sixteen
	};

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	SersicProfile(const Model &model, const std::string &name);

	void validate() override;

protected:

	/*
	 * ----------------------
	 * Inherited from Profile
	 * ----------------------
	 */
	bool parameter_impl(const std::string &name, double val) override;
	bool parameter_impl(const std::string &name, bool val) override;

	/*
	 * ----------------------------
	 * Inherited from RadialProfile
	 * ----------------------------
	 */
	void evaluate(std::vector<double> &image) override;
	void initial_calculations() override;
	void subsampling_params(double x, double y, unsigned int &res, unsigned int &max_rec) override;
	double get_pixel_scale() override;

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
	 * The effective radius
	 */
	double re;

	/**
	 * The sersic index
	 */
	double nser;

	/**
	 * Rescale flux up to rscale_max or not
	 */
	bool rescale_flux;
	// @}

	/* these are internally calculated when the profile is evaluated */
	double _bn;
	double _rescale_factor;

private:

	double (*m_eval_function)(double x, double y, double box, double re, double nser, double bn);

	template <bool boxy, SersicProfile::rfactor_invexp_t t>
	void init_eval_function();

	double fluxfrac(double fraction) const;

#ifdef PROFIT_OPENCL

protected:
	virtual void add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const override;
	virtual void add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const override;
	virtual bool supports_opencl() const override;

private:
	template <typename FT>
	void add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const;

#endif /* PROFIT_OPENCL */

};

} /* namespace profit */

#endif /* PROFIT_SERSIC_H */
