/**
 * CoreSersic profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Rodrigo Tobar
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
#include "profit/coresersic.h"
#include "profit/exceptions.h"
#include "profit/utils.h"


namespace profit
{

/**
 * The evaluation of the coresersic profile at coresersic coordinates (x,y).
 *
 * The coresersic profile has this form:
 *
 *    (1+(r/rb)^(-a))^(b/a)*
 *        exp(-bn*(((r^a+rb^a)/re^a))^(1/(nser*a)))
 *
 * with::
 *
 *           r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *           B = box parameter
 */
double CoreSersicProfile::evaluate_at(double x, double y) const {

	using std::abs;
	using std::exp;
	using std::pow;

	double box = this->box + 2.;
	double r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	return pow(1 + pow(r/rb,-a), b/a) *
	       exp(-_bn * pow((pow(r, a) + pow(rb, a))/pow(re,a), 1/(nser*a)));
}

void CoreSersicProfile::validate() {

	RadialProfile::validate();

	if ( re <= 0 ) {
		throw invalid_parameter("re <= 0, must have re > 0");
	}
	if ( rb <= 0 ) {
		throw invalid_parameter("rb <= 0, must have rb > 0");
	}
	if ( nser <= 0 ) {
		throw invalid_parameter("nser <= 0, must have nser > 0");
	}
	if ( a <= 0 ) {
		throw invalid_parameter("a <= 0, must have a > 0");
	}
	if ( b > 1.999 ) {
		throw invalid_parameter("b > 1.999, must have b < 1.999");
	}

}

double CoreSersicProfile::integrate_at(double r) const {

	using std::exp;
	using std::pow;

	return r * pow(1 + pow(r/rb,-a), b/a) *
	       exp(-_bn * pow((pow(r, a) + pow(rb, a))/pow(re,a), 1/(nser*a)));
}

double CoreSersicProfile::get_lumtot(double r_box) {

	/*
	 * We numerically integrate r from 0 to infinity
	 * to get the total luminosity
	 */
	auto int_f = [](double r, void *ctx){
		CoreSersicProfile *p = static_cast<CoreSersicProfile *>(ctx);
		return p->integrate_at(r);
	};
	double magtot = integrate_qagi(int_f, 0, this);
	return 2 * M_PI * axrat * magtot/r_box;
}

void CoreSersicProfile::initial_calculations() {

	/*
	 * bn needs to be calculated before calling the super method
	 * because it's used to calculate the total luminosity
	 */
	this->_bn = qgamma(0.5, 2*this->nser);

	/* Common calculations first */
	RadialProfile::initial_calculations();

}

double CoreSersicProfile::adjust_rscale_switch() {
	return 1;
}

double CoreSersicProfile::adjust_rscale_max() {
	return 100;
}

double CoreSersicProfile::get_rscale() {
	return this->re;
}

double CoreSersicProfile::adjust_acc() {
	return this->acc;
}

CoreSersicProfile::CoreSersicProfile(const Model &model, const std::string &name) :
	RadialProfile(model, name),
	re(1), rb(1), nser(4), a(1), b(1)
{
	// no-op
}
bool CoreSersicProfile::parameter_impl(const std::string &name, double val) {

	if( RadialProfile::parameter_impl(name, val) ) {
		return true;
	}

	if( name == "re" )        { re = val; }
	else if( name == "rb" )   { rb = val; }
	else if( name == "nser" ) { nser = val; }
	else if( name == "a" )    { a = val; }
	else if( name == "b" )    { b = val; }
	else {
		return false;
	}

	return true;
}

#ifdef PROFIT_OPENCL
void CoreSersicProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<float>(index, kernel);
}

void CoreSersicProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<double>(index, kernel);
}

template <typename FT>
void CoreSersicProfile::add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const {
	kernel.setArg(index++, static_cast<FT>(re));
	kernel.setArg(index++, static_cast<FT>(rb));
	kernel.setArg(index++, static_cast<FT>(nser));
	kernel.setArg(index++, static_cast<FT>(a));
	kernel.setArg(index++, static_cast<FT>(b));
	kernel.setArg(index++, static_cast<FT>(_bn));
}

#endif /* PROFIT_OPENCL */

} /* namespace profit */