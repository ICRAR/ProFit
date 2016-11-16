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

#include "profit/brokenexponential.h"
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
static
double _brokenexponential_for_xy_r(const RadialProfile &sp,
                                   double x, double y,
                                   double r, bool reuse_r) {

	const BrokenExponentialProfile &bep = static_cast<const BrokenExponentialProfile &>(sp);
	if( !reuse_r && bep.box == 0 ) {
		r = sqrt(x*x + y*y);
	}
	else if( !reuse_r ) { // && bep.box != 0
		double box = bep.box + 2.;
		r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}
	// else reuse_r == true, so r is used as is

	return exp(-r/bep.h1)*pow(1+exp(bep.a*(r-bep.rb)),(1/bep.a)*(1/bep.h1-1/bep.h2));

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

eval_function_t BrokenExponentialProfile::get_evaluation_function() {
	return &_brokenexponential_for_xy_r;
}

static
double brokenexponential_int(double r, void *ex) {
	BrokenExponentialProfile *bep = (BrokenExponentialProfile *)ex;
	return r * exp(-r/bep->h1)*pow(1+exp(bep->a*(r-bep->rb)),(1/bep->a)*(1/bep->h1-1/bep->h2));
}

double BrokenExponentialProfile::get_lumtot(double r_box) {
	/*
	 * We numerically integrate r from 0 to infinity
	 * to get the total luminosity
	 */
	double magtot = integrate_qagi(&brokenexponential_int, 0, this);
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

BrokenExponentialProfile::BrokenExponentialProfile(const Model &model) :
	RadialProfile(model),
	h1(1), h2(1), rb(1), a(1)
{
	// no-op
}

} /* namespace profit */
