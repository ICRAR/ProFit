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

#include "profit/coresersic.h"
#include "profit/utils.h"

using namespace std;

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
static
double _coresersic_for_xy_r(const RadialProfile &sp,
                            double x, double y,
                            double r, bool reuse_r) {

	const CoreSersicProfile &csp = static_cast<const CoreSersicProfile &>(sp);

	if( csp.box == 0 && !reuse_r ) {
		r = sqrt(x*x + y*y);
	}
	else if( csp.box != 0 ){
		double box = csp.box + 2.;
		r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}
	// else csp.box == 0 && reuse_r, so we leave r untouched

	double rb = csp.rb;
	double a = csp.a;
	double b = csp.b;
	double bn = csp._bn;
	double re = csp.re;
	double nser = csp.nser;

	return pow(1 + pow(r/rb,-a), b/a) *
	       exp(-bn * pow((pow(r, a) + pow(rb, a))/pow(re,a), 1/(nser*a)));

}

eval_function_t CoreSersicProfile::get_evaluation_function() {
	return &_coresersic_for_xy_r;
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

static
double coresersic_int(double r, void *ex) {
	CoreSersicProfile *csp = (CoreSersicProfile *)ex;
	double rb = csp->rb;
	double a = csp->a;
	double b = csp->b;
	double bn = csp->_bn;
	double re = csp->re;
	double nser = csp->nser;
	return r * pow(1 + pow(r/rb,-a), b/a) *
	       exp(-bn * pow((pow(r, a) + pow(rb, a))/pow(re,a), 1/(nser*a)));
}


double CoreSersicProfile::get_lumtot(double r_box) {

	/*
	 * We numerically integrate r from 0 to infinity
	 * to get the total luminosity
	 */
	double magtot = integrate_qagi(&coresersic_int, 0, this);
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

CoreSersicProfile::CoreSersicProfile(const Model &model) :
	RadialProfile(model),
	re(1), rb(1), nser(4), a(1), b(1)
{
	// no-op
}

} /* namespace profit */
