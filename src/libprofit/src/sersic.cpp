/**
 * Sersic profile implementation
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
#include "profit/model.h"
#include "profit/opencl.h"
#include "profit/sersic.h"
#include "profit/utils.h"

namespace profit
{

/*
 * The evaluation of the sersic profile at sersic coordinates (x,y).
 *
 * The sersic profile has this form:
 *
 * e^{-bn * (r_factor - 1)}
 * where r_factor = (r/Re)^{1/nser}
 *              r = (x^{2+b} + y^{2+b})^{1/(2+b)}
 *              b = box parameter
 *
 * Reducing:
 *  r_factor = ((x/re)^{2+b} + (y/re)^{2+b})^{1/(nser*(2+b)}
 *
 * Although this reduced form has only three powers instead of 4
 * in the original form, it still doesn't mean that it's the fastest
 * way of computing r_factor. Depending on the libc/libm being used
 * using a combination of sqrt, cbrt and pow will yield better or
 * worse performances (within certain limits).
 *
 * Also, in particular for b = 0:
 *  r_factor = ((x*x + y*y/(re*re))^{1/(2*nser}
 *
 * In our code we thus logically decompose r_factor (for both cases) as such:
 *
 * r_factor = base^(1/invexp)
 *
 * We then use specialized code templates to get the proper base, the proper
 * invexp, and finally to check whether some calls to pow() can be avoided or
 * not.
 */

/*
 * r_factor calculation follows. Several template specializations avoid the
 * call to pow().
 * This first generic template will finally be called only with parameters
 * true/general because all the other combinations are already specialized.
 */

/* "true" cases for r_factor */
template<bool boxy, SersicProfile::rfactor_invexp_t t>
inline double _r_factor(double b, double invexp)
{
  return std::pow(b, 1/invexp);
}

template<> inline double _r_factor<true, SersicProfile::pointfive>(double b, double  /*invexp*/)
{
	return b*b;
}

template<> inline double _r_factor<true, SersicProfile::one>(double b, double  /*invexp*/)
{
	return b;
}

template<> inline double _r_factor<true, SersicProfile::two>(double b, double  /*invexp*/)
{
	return std::sqrt(b);
}

template<> inline double _r_factor<true, SersicProfile::three>(double b, double  /*invexp*/)
{
	return std::cbrt(b);
}

template<> inline double _r_factor<true, SersicProfile::four>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(b));
}

template<> inline double _r_factor<true, SersicProfile::eight>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(std::sqrt(b)));
}

template<> inline double _r_factor<true, SersicProfile::sixteen>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(std::sqrt(std::sqrt(b))));
}

/* "false" cases for r_factor */
template<> inline double _r_factor<false, SersicProfile::general>(double b, double invexp)
{
	return std::pow(std::sqrt(b), 1/invexp);
}

template<> inline double _r_factor<false, SersicProfile::pointfive>(double b, double  /*invexp*/)
{
	return b;
}

template<> inline double _r_factor<false, SersicProfile::one>(double b, double  /*invexp*/)
{
	return std::sqrt(b);
}

template<> inline double _r_factor<false, SersicProfile::two>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(b));
}

template<> inline double _r_factor<false, SersicProfile::three>(double b, double  /*invexp*/)
{
	return std::cbrt(std::sqrt(b));
}

template<> inline double _r_factor<false, SersicProfile::four>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(std::sqrt(b)));
}

template<> inline double _r_factor<false, SersicProfile::eight>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(std::sqrt(std::sqrt(b))));
}

template<> inline double _r_factor<false, SersicProfile::sixteen>(double b, double  /*invexp*/)
{
	return std::sqrt(std::sqrt(std::sqrt(std::sqrt(std::sqrt(b)))));
}

/*
 * The base component of r_factor
 */
template<bool boxy>
inline double _base(double x, double y, double re, double exponent)
{
	return std::pow(std::abs(x/re), exponent) + std::pow(std::abs(y/re), exponent);
}

template<>
inline double _base<false>(double x, double y, double re, double  /*exponent*/)
{
	return (x*x + y*y)/(re * re);
}

/*
 * The invexpt component of r_factor
 */
template<bool boxy>
inline double _invexp(const double nser, const double exponent)
{
  return nser*exponent;
}

template<>
inline double _invexp<false>(const double nser, const double  /*exponent*/)
{
  return nser;
}

template<bool boxy, SersicProfile::rfactor_invexp_t t>
static inline
double eval_function(double x, double y, double box, double re, double nser, double bn) {
	double B = box + 2;
	double base = _base<boxy>(x, y, re, B);
	double r_factor = _r_factor<boxy, t>(base, _invexp<boxy>(nser, B));
	return std::exp(-bn * (r_factor - 1));
}

/*
 * The main sersic evaluation function for a given X/Y coordinate
 */
double SersicProfile::evaluate_at(double x, double y) const {
	return m_eval_function(x, y, box, re, nser, _bn);
}

void SersicProfile::validate() {

	RadialProfile::validate();

	if ( re <= 0 ) {
		throw invalid_parameter("re <= 0, must have re > 0");
	}
	if ( nser <= 0 ) {
		throw invalid_parameter("nser <= 0, must have nser > 0");
	}

}

template <bool boxy, SersicProfile::rfactor_invexp_t t>
void SersicProfile::init_eval_function() {
	m_eval_function = eval_function<boxy, t>;
}

void SersicProfile::evaluate(Image &image, const Mask &mask, const PixelScale &scale,
    const Point &offset, double magzero) {

	// inv_exponent is exactly what is yield by the templated _invexp function
	// later on during each individual evaluation
	// We need to calculate it here though to decide which template to choose
	double inv_exponent = nser;
	if( this->box != 0 ) {
		inv_exponent *= (box + 2);
		if( almost_equals(inv_exponent, 0.5) )    init_eval_function<true, pointfive>();
		else if( almost_equals(inv_exponent, 1) ) init_eval_function<true, one>();
		else if( almost_equals(inv_exponent, 2) ) init_eval_function<true, two>();
		else if( almost_equals(inv_exponent, 3) ) init_eval_function<true, three>();
		else if( almost_equals(inv_exponent, 4) ) init_eval_function<true, four>();
		else if( almost_equals(inv_exponent, 8) ) init_eval_function<true, eight>();
		else if( almost_equals(inv_exponent, 16)) init_eval_function<true, sixteen>();
		else                                      init_eval_function<true, general>();
	}
	else {
		if( almost_equals(inv_exponent, 0.5) )     init_eval_function<false, pointfive>();
		else if( almost_equals(inv_exponent, 1) )  init_eval_function<false, one>();
		else if( almost_equals(inv_exponent, 2) )  init_eval_function<false, two>();
		else if( almost_equals(inv_exponent, 3) )  init_eval_function<false, three>();
		else if( almost_equals(inv_exponent, 4) )  init_eval_function<false, four>();
		else if( almost_equals(inv_exponent, 8) )  init_eval_function<false, eight>();
		else if( almost_equals(inv_exponent, 16) ) init_eval_function<false, sixteen>();
		else                                       init_eval_function<false, general>();
	}

	return RadialProfile::evaluate(image, mask, scale, offset, magzero);
}

double SersicProfile::fluxfrac(double fraction) const {
	double ratio = qgamma(fraction, 2*nser) / _bn;
	return re * std::pow(ratio, nser);
}

double SersicProfile::adjust_rscale_switch() {
	/*
	 * Find the point at which we capture most of the flux (sensible place
	 * for upscaling). We make sure upscaling doesn't go beyond 20 pixels,
	 * but don't let it become less than 1 pixel (means we do no worse than
	 * GALFIT anywhere)
	 */
	double rscale_switch = std::ceil(fluxfrac(1. - nser*nser/2e3));
	rscale_switch = std::max(std::min(rscale_switch, 20.), 2.);
	return rscale_switch / this->re;
}

double SersicProfile::adjust_rscale_max() {
	return std::ceil(std::max(fluxfrac(0.9999), 2.) / re);
}

double SersicProfile::adjust_acc(double acc) {
	acc /= std::sqrt(nser);
	return std::max(0.1, acc) / this->axrat;
}

double SersicProfile::get_lumtot() {
	using std::exp;
	using std::pow;
	double g_factor = gammafn(2*this->nser);
	return pow(this->re, 2) * 2 * M_PI * this->nser * g_factor *
	       exp(this->_bn)/pow(this->_bn, 2*this->nser);
}

void SersicProfile::initial_calculations() {

	/*
	 * bn needs to be calculated before calling the super method
	 * because it's used to calculate the total luminosity
	 */
	this->_bn = qgamma(0.5, 2*this->nser);

	/* Common calculations first */
	RadialProfile::initial_calculations();

	/* Just some additional adjustments on rescale_factor */
	if( this->adjust ) {
		this->_rescale_factor = 1;
		if( this->rescale_flux ) {
			double flux_r;
			flux_r = this->_bn * std::pow(this->rscale_max/this->re, 1/this->nser);
			flux_r = pgamma(flux_r, 2*this->nser);
			this->_rescale_factor = 1/flux_r;
		}
	}

}

/**
 * The scale by which each image pixel value is multiplied.
 * The sersic profile supports a rescale factor that is applied here.
 */
double SersicProfile::get_pixel_scale(const PixelScale &scale) {
	double flux_scale = RadialProfile::get_pixel_scale(scale);
	if( this->rescale_flux ) {
		flux_scale *= this->_rescale_factor;
	}
	return flux_scale;
}

double SersicProfile::get_rscale() {
	return this->re;
}

void SersicProfile::subsampling_params(double x, double y,
                                       unsigned int &resolution,
                                       unsigned int &max_recursions) {

	using std::abs;

	RadialProfile::subsampling_params(x, y, resolution, max_recursions);

	/* Higher subsampling params for central pixel if nser > 1 (only when auto-adjusting) */
	bool center_pixel = abs(x - _xcen) < model.get_image_pixel_scale().first &&
	                    abs(y - _ycen) < model.get_image_pixel_scale().second;
	if( adjust && center_pixel && nser > 1 ) {
		resolution = 8;
		max_recursions = 10;
	}

}

/**
 * The sersic creation function
 */
SersicProfile::SersicProfile(const Model &model, const std::string &name) :
	RadialProfile(model, name),
	re(1), nser(1),
	rescale_flux(false)
{
	register_parameter("re", re);
	register_parameter("nser", nser);
	register_parameter("rescale_flux", rescale_flux);
}

#ifdef PROFIT_OPENCL
void SersicProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<float>(index, kernel);
}

void SersicProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	add_kernel_parameters<double>(index, kernel);
}

template <typename FT>
void SersicProfile::add_kernel_parameters(unsigned int index, cl::Kernel &kernel) const {
	kernel.setArg((index++), FT(nser));
	kernel.setArg((index++), FT(_bn));
}

#endif /* PROFIT_OPENCL */

} /* namespace profit */
