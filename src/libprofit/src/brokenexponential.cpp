/**
 * BrokenExponential profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham
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

#include <cmath>

#include "profit/common.h"
#include "profit/brokenexponential.h"
#include "profit/exceptions.h"
#include "profit/utils.h"


namespace profit
{

inline static
double _broken_exponential(double r, double h1, double h2, double rb, double a) {

	/*
	 * The broken exponential profile for radius r is:
	 *
	 *  exp(-r/h1)*(1+exp(a*(r-rb)))^((1/a)*(1/h1-1/h2))
	 *
	 * The problem with this direct approach is that exp(base) diverges whereas
	 * exp(-r/h1) converges to zero. To avoid this we perform the following
	 * replacements and rewrite the equation:
	 *
	 *   base = r - rb
	 *   exponent = 1/h1 - 1/h2
	 *
	 *    exp(-r/h1) * (1 + exp(a * base)) ^ (exponent / a)
	 *  = exp(-r/h1) * exp(log((1 + exp(a * base)) ^ (exponent / a)))
	 *  = exp(-r/h1 + log((1 + exp(a * base)) ^ (exponent / a))
	 *  = exp(-r/h1 + exponent / a * log(1 + exp(a * base)))
	 *
	 * In this last expression, when (a * base) is big, then doing
	 * log(1 + exp(a * base)) yields the same result as log(exp(a * base))
	 * (which equals a * base, of course). This happens already at a * base = 34,
	 * although we check 40 just to be conservative.
	 * Thus, the final result in this case becomes:
	 *
	 *  = exp(-r/h1 + exponent * base)
	 */

	using std::exp;
	using std::log;

	double base = r - rb;
	double expo = 1 / h1 - 1 / h2;
	if (a * base < 40) {
		base = log(1 + exp(a * base)) / a;
	}

	return exp(-r / h1 + expo * base);
}

/**
 * The evaluation of the brokenexponential profile at brokenexponential coordinates (x,y).
 *
 * The broken exponential profile has this form:
 *
 *    exp(-r/h1)*(1+exp(a*(r-rb)))^((1/a)*(1/h1-1/h2))
 *
 * with::
 *
 *       r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *       B = box parameter
 */
double BrokenExponentialProfile::evaluate_at(double x, double y) const {

	using std::abs;
	using std::pow;

	double box = this->box + 2.;
	double r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	return _broken_exponential(r, h1, h2, rb, a);
}

void BrokenExponentialProfile::validate() {

	RadialProfile::validate();

	if ( h1 <= 0 ) {
		throw invalid_parameter("h1 <= 0, must have h1 > 0");
	}
	if ( h2 <= 0 ) {
		throw invalid_parameter("h2 <= 0, must have h2 > 0");
	}
	if ( rb <= 0 ) {
		throw invalid_parameter("rb <= 0, must have rb > 0");
	}
}

double BrokenExponentialProfile::integrate_at(double r) const {
	return r * _broken_exponential(r, h1, h2, rb, a);
}

double BrokenExponentialProfile::get_lumtot(double r_box) {
	/*
	 * We numerically integrate r from 0 to infinity
	 * to get the total luminosity
	 */
	double magtot = integrate_qagi(
		[](double r, void *ctx) -> double {
			BrokenExponentialProfile *p = static_cast<BrokenExponentialProfile *>(ctx);
			return p->integrate_at(r);
		},
		0, this);
	return 2*M_PI * axrat * magtot/r_box;
}

double BrokenExponentialProfile::adjust_rscale_switch() {
	return 1;
}

double BrokenExponentialProfile::adjust_rscale_max() {
	return 100;
}

double BrokenExponentialProfile::get_rscale() {
	return this->h1*4;
}

double BrokenExponentialProfile::adjust_acc() {
	return this->acc;
}

BrokenExponentialProfile::BrokenExponentialProfile(const Model &model, const std::string &name) :
	RadialProfile(model, name),
	h1(1), h2(1), rb(1), a(1)
{
	// no-op
}

bool BrokenExponentialProfile::parameter_impl(const std::string &name, double val) {

	if( RadialProfile::parameter_impl(name, val) ) {
		return true;
	}

	if( name == "h1" )      { h1 = val; }
	else if( name == "h2" ) { h2 = val; }
	else if( name == "rb" ) { rb = val; }
	else if( name == "a" )  { a = val; }
	else {
		return false;
	}

	return true;
}

#ifdef PROFIT_OPENCL
void BrokenExponentialProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<float>(index, kernel);
}

void BrokenExponentialProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<double>(index, kernel);
}

template <typename FT>
void BrokenExponentialProfile::add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const {
	kernel.setArg(index++, static_cast<FT>(h1));
	kernel.setArg(index++, static_cast<FT>(h2));
	kernel.setArg(index++, static_cast<FT>(rb));
	kernel.setArg(index++, static_cast<FT>(a));
}

#endif /* PROFIT_OPENCL */

} /* namespace profit */