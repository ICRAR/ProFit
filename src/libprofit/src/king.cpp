/**
 * King profile implementation
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
#include "profit/exceptions.h"
#include "profit/king.h"
#include "profit/utils.h"


namespace profit
{

/*
 * The evaluation of the king profile at king coordinates (x,y).
 *
 * The king profile has this form:
 *
 *    (1-temp)^(-a) * (1/(1+(r/rc)^2)^(1/a)-temp)^a
 *
 * with::
 *
 *   temp  = 1/(1+(rt/rc)^2)^(1/a)
 *       r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *       B = box parameter
 */
double KingProfile::evaluate_at(double x, double y) const {

	using std::abs;
	using std::pow;

	double box = this->box + 2.;
	double r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);

	if( r < rt ) {
		return pow(1/pow(1 + pow(r/rc, 2), 1/a) - 1/pow(1 + pow(rt/rc, 2), 1/a), a);
	}

	return 0;
}

void KingProfile::validate() {

	RadialProfile::validate();

	if ( rc <= 0 ) {
		throw invalid_parameter("rc <= 0, must have rc > 0");
	}
	if ( rt <= 0 ) {
		throw invalid_parameter("rt <= 0, must have rt > 0");
	}
	if ( a < 0 ) {
		throw invalid_parameter("a < 0, must have a >=0");
	}

}

double KingProfile::integrate_at(double r) const {

	using std::pow;

	if( r < rt ) {
		return r * pow(1/pow(1 + pow(r/rc, 2), 1/a) - 1/pow(1 + pow(rt/rc, 2), 1/a), a);
	}
	return 0;
}

double KingProfile::get_lumtot(double r_box) {
	/*
	 * We numerically integrate r from 0 to rt
	 * to get the total luminosity
	 */
	double magtot = integrate_qags([](double r, void *ctx) -> double {
		auto kp = static_cast<KingProfile *>(ctx);
		return kp->integrate_at(r);
	}, 0, this->rt, this);
	return 2*M_PI * axrat * magtot/r_box;
}

double KingProfile::adjust_rscale_switch() {
	return 0.5;
}

double KingProfile::adjust_rscale_max() {
	return 1.5;
}

double KingProfile::get_rscale() {
	return this->rt;
}

double KingProfile::adjust_acc() {
	return this->acc;
}

KingProfile::KingProfile(const Model &model, const std::string &name) :
	RadialProfile(model, name),
	rc(1), rt(3), a(2)
{
	// no-op
}

bool KingProfile::parameter_impl(const std::string &name, double val) {

	if( RadialProfile::parameter_impl(name, val) ) {
		return true;
	}

	if( name == "rc" )      { rc = val; }
	else if( name == "rt" ) { rt = val; }
	else if( name == "a" )  { a = val; }
	else {
		return false;
	}

	return true;
}

#ifdef PROFIT_OPENCL
void KingProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<float>(index, kernel);
}

void KingProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<double>(index, kernel);
}

template <typename FT>
void KingProfile::add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const {
	kernel.setArg(index++, static_cast<FT>(rc));
	kernel.setArg(index++, static_cast<FT>(rt));
	kernel.setArg(index++, static_cast<FT>(a));
}

#endif /* PROFIT_OPENCL */

} /* namespace profit */