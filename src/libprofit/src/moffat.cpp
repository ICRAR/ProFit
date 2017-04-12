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

#include "profit/common.h"
#include "profit/exceptions.h"
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
double MoffatProfile::evaluate_at(double x, double y) const {
	double box = 2 + this->box;
	double r = pow( pow(abs(x), box) + pow(abs(y), box), 1./(box));
	double r_factor = r/rscale;
	return pow(1 + r_factor*r_factor, -con);
}

void MoffatProfile::validate() {

	RadialProfile::validate();

	if ( fwhm <= 0 ) {
		throw invalid_parameter("fwhm <= 0, must have fwhm > 0");
	}
	if ( con < 0 ) {
		throw invalid_parameter("con < 0, must have con >= 0");
	}

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
	return 8;
}

double MoffatProfile::adjust_acc() {
	return 0.1/axrat;
}


MoffatProfile::MoffatProfile(const Model &model, const string &name) :
	RadialProfile(model, name),
	fwhm(3), con(2)
{
	// no-op
}

bool MoffatProfile::parameter_impl(const string &name, double val) {

	if( RadialProfile::parameter_impl(name, val) ) {
		return true;
	}

	if( name == "fwhm" )     { fwhm = val; }
	else if( name == "con" ) { con = val; }
	else {
		return false;
	}

	return true;
}

#ifdef PROFIT_OPENCL
void MoffatProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<float>(index, kernel);
}

void MoffatProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<double>(index, kernel);
}

template <typename FT>
void MoffatProfile::add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const {
	kernel.setArg(index++, static_cast<FT>(con));
}

#endif /* PROFIT_OPENCL */

} /* namespace profit */
