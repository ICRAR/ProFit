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


using namespace std;

namespace profit
{

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
	double box = this->box + 2.;
	double r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	return exp(-r/h1)*pow(1+exp(a*(r-rb)),(1/a)*(1/h1-1/h2));

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
	if ( h2 > h1 ) {
		throw invalid_parameter("h2 > h1, must have h2 <= h1");
	}

}

double BrokenExponentialProfile::integrate_at(double r) const {
	return r * exp(-r/h1)*pow(1+exp(a*(r-rb)),(1/a)*(1/h1-1/h2));
}

double BrokenExponentialProfile::get_lumtot(double r_box) {
	/*
	 * We numerically integrate r from 0 to infinity
	 * to get the total luminosity
	 */
	double magtot = integrate_qagi(
		[](double r, void *ctx) -> double {
			BrokenExponentialProfile *p = reinterpret_cast<BrokenExponentialProfile *>(ctx);
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

BrokenExponentialProfile::BrokenExponentialProfile(const Model &model, const string &name) :
	RadialProfile(model, name),
	h1(1), h2(1), rb(1), a(1)
{
	// no-op
}

bool BrokenExponentialProfile::parameter_impl(const string &name, double val) {

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

} /* namespace profit */