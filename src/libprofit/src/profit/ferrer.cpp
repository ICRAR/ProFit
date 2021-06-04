/**
 * Ferrer profile implementation
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
#include "profit/ferrer.h"
#include "profit/utils.h"


namespace profit
{

/**
 * The evaluation of the ferrer profile at ferrer coordinates (x,y).
 *
 * The ferrer profile has this form:
 *
 * (1-r_factor)^(a)
 *
 * where r_factor = (r/rscale)^(2-b)
 *              r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *              B = box parameter
 */
double FerrerProfile::evaluate_at(double x, double y) const {

	using std::pow;

	double r = boxy_r(x, y);
	double r_factor = r/rscale;
	if( r_factor < 1 ) {
		return pow(1 - pow(r_factor, 2 - b), a);
	}
	return 0;
}

void FerrerProfile::validate() {

	RadialProfile::validate();

	if ( rout <= 0 ) {
		throw invalid_parameter("rout <= 0, must have rout >= 0");
	}
	if ( a < 0 ) {
		throw invalid_parameter("a < 0, must have a >= 0");
	}
	if ( b > 2 ) {
		throw invalid_parameter("b > 2, must have b <= 2");
	}

}

double FerrerProfile::get_lumtot() {

	using std::pow;

	/*
	 * Wolfram Alpha gave for g_factor:
	 *
	 *   G(a + 1) * G((4-b)/(2-b)) / G(a + 2/(2-b) + 1)
	 *
	 * But (4-b)/(2-b) == 1 + 2/(2-b), and G(a+1) == a*G(a)
	 * Thus the expression above equals to:
	 *
	 *   a * G(a) * G(1 + 2/(2-b)) / G(a + 2/(2-b) + 1)
	 *   a * B(a, 1+2/(2-b))
	 *
	 * Although the results are the same, high-level mathematical libs
	 * like GSL and R handle calculation errors better than we do and
	 * thus it's better to let them deal with intermediate results.
	 * For example b=1.99 gives NaN with our gamma-based calculations,
	 * but still converges using beta.
	 */
	double g_factor = a * beta(a, 1 + 2/(2-b));
	return pow(rout, 2) * M_PI * g_factor;
}

double FerrerProfile::adjust_rscale_switch() {
	return 0.5;
}

double FerrerProfile::adjust_rscale_max() {
	return 1;
}

double FerrerProfile::get_rscale() {
	return this->rout;
}

FerrerProfile::FerrerProfile(const Model &model, const std::string &name) :
	RadialProfile(model, name),
	rout(3), a(1), b(1)
{
	// this profile defaults to a different accuracy
	this->acc = 1;

	register_parameter("rout", rout);
	register_parameter("a", a);
	register_parameter("b", b);
}

#ifdef PROFIT_OPENCL
void FerrerProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<float>(index, kernel);
}

void FerrerProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<double>(index, kernel);
}

template <typename FT>
void FerrerProfile::add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const {
	kernel.setArg((index++), FT(a));
	kernel.setArg((index++), FT(b));
}

#endif /* PROFIT_OPENCL */

} /* namespace profit */
