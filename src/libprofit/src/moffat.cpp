/**
 * Moffat profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Dan Taranu, Rodrigo Tobar
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
#include <algorithm>

#include "profit/moffat.h"

using namespace std;

namespace profit
{

/*
 * The evaluation of the moffat profile at moffat coordinates (x,y).
 *
 * The moffat profile has this form:
 *
 * (1+r_factor^2)^(-c)
 *
 * where r_factor = (r/rscale)
 *              r = (x^{2+b} + y^{2+b})^{1/(2+b)}
 *              b = box parameter
 *
 * Reducing:
 *  r_factor = ((x/rscale)^{2+b} + (y/rscale)^{2+b})^{1/(2+b)}
 */
static
double _moffat_for_xy_r(RadialProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	MoffatProfile *mp = static_cast<MoffatProfile *>(sp);
	double r_factor;
	if( mp->box == 0 ) {
		r_factor = sqrt(x*x + y*y);
	}
	else {
		double box = 2 + mp->box;
		r_factor = pow( pow(abs(x), box) + pow(abs(y), box), 1./(box));
	}

	r_factor /= mp->rscale;
	return pow(1 + r_factor*r_factor, -mp->con);
}

eval_function_t MoffatProfile::get_evaluation_function() {
	return &_moffat_for_xy_r;
}

double MoffatProfile::get_lumtot(double r_box) {
	double con = this->con;
	return pow(this->rscale, 2) * M_PI * axrat/(con-1)/r_box;
}

double MoffatProfile::get_rscale() {
	return fwhm/(2*sqrt(pow(2,(1/con))-1));
}

double MoffatProfile::adjust_rscale_switch() {
	double rscale_switch = this->fwhm*4;
	rscale_switch = max(min(rscale_switch, 20.), 2.);
	return rscale_switch / this->rscale;
}

double MoffatProfile::adjust_rscale_max() {
	return this->fwhm*8;
}

double MoffatProfile::adjust_acc() {
	return 0.1/axrat;
}


MoffatProfile::MoffatProfile() :
	RadialProfile(),
	fwhm(3), con(2)
{
	// no-op
}

} /* namespace profit */
