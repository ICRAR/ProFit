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

#include "profit/king.h"
#include "profit/utils.h"

using namespace std;

namespace profit
{

/**
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
static
double _king_for_xy_r(const RadialProfile &sp,
                      double x, double y,
                      double r, bool reuse_r) {

	const KingProfile &kp = static_cast<const KingProfile &>(sp);
	if( !reuse_r && kp.box == 0 ) {
		r = sqrt(x*x + y*y);
	}
	else if( !reuse_r ) { // && kp.box != 0
		double box = kp.box + 2.;
		r = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}
	// else reuse_r == true, so r is used as is

	if( r < kp.rt ) {
		return pow(1/pow(1 + pow(r/kp.rc, 2), 1/kp.a) - 1/pow(1 + pow(kp.rt/kp.rc, 2), 1/kp.a), kp.a);
	}

	return 0;
}

void KingProfile::validate() {

	RadialProfile::validate();

	if ( rc <= 0 ) {
		throw invalid_parameter("rc <= 0, must have rc > 0");
	}
	if ( rt <= 0 ) {
		throw invalid_parameter("rt <= 0, must have rc > 0");
	}
	if ( a < 0 ) {
		throw invalid_parameter("a < 0, must have a >=0");
	}

}

eval_function_t KingProfile::get_evaluation_function() {
	return &_king_for_xy_r;
}

static
double king_int(double r, void *ex) {
	KingProfile *kp = (KingProfile *)ex;
	if( r < kp->rt ) {
		return r * pow(1/pow(1 + pow(r/kp->rc, 2), 1/kp->a) - 1/pow(1 + pow(kp->rt/kp->rc, 2), 1/kp->a), kp->a);
	}
	return 0;
}

double KingProfile::get_lumtot(double r_box) {
	/*
	 * We numerically integrate r from 0 to rt
	 * to get the total luminosity
	 */
	double magtot = integrate_qags(&king_int, 0, this->rt, this);
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

KingProfile::KingProfile(const Model &model) :
	RadialProfile(model),
	rc(1), rt(3), a(2)
{
	// no-op
}

} /* namespace profit */
