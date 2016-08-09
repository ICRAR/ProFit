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

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

/*
 * We use either GSL or Rmath to provide the low-level
 * beta, gamma and qgamma_inv functions needed by the sersic profile.
 * If neither is given, the user will have to feed the profiles with
 * the appropriate function pointers after creating them.
 */
#if defined(HAVE_GSL)
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#elif defined(HAVE_R)
#define R_NO_REMAP
#include <Rmath.h>
#endif

#include "sersic.h"

using namespace std;

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
 * The nser parameter is a double; we need an enumeration of the known values
 * to optimize for to use in our templates
 */
enum nser_t {
	general,
	pointfive,
	one,
	two,
	three,
	four,
	eight,
	sixteen
};

/*
 * r_factor calculation follows. Several template specializations avoid the
 * call to pow().
 * This first generic template will finally be called only with parameters
 * true/general because all the other combinations are already specialized.
 */

/* "true" cases for r_factor */
template<bool boxy, nser_t t>
inline double _r_factor(double b, double invexp)
{
  return pow(b, 1/invexp);
}

template<> inline double _r_factor<true, pointfive>(double b, double invexp)
{
	return b*b;
}

template<> inline double _r_factor<true, one>(double b, double invexp)
{
	return b;
}

template<> inline double _r_factor<true, two>(double b, double invexp)
{
	return sqrt(b);
}

template<> inline double _r_factor<true, three>(double b, double invexp)
{
	return cbrt(b);
}

template<> inline double _r_factor<true, four>(double b, double invexp)
{
	return sqrt(sqrt(b));
}

template<> inline double _r_factor<true, eight>(double b, double invexp)
{
	return sqrt(sqrt(sqrt(b)));
}

template<> inline double _r_factor<true, sixteen>(double b, double invexp)
{
	return sqrt(sqrt(sqrt(sqrt(b))));
}

/* "false" cases for r_factor */
template<> inline double _r_factor<false, general>(double b, double invexp)
{
	return pow(sqrt(b), 1/invexp);
}

template<> inline double _r_factor<false, pointfive>(double b, double invexp)
{
	return b;
}

template<> inline double _r_factor<false, one>(double b, double invexp)
{
	return sqrt(b);
}

template<> inline double _r_factor<false, two>(double b, double invexp)
{
	return sqrt(sqrt(b));
}

template<> inline double _r_factor<false, three>(double b, double invexp)
{
	return cbrt(sqrt(b));
}

template<> inline double _r_factor<false, four>(double b, double invexp)
{
	return sqrt(sqrt(sqrt(b)));
}

template<> inline double _r_factor<false, eight>(double b, double invexp)
{
	return sqrt(sqrt(sqrt(sqrt(b))));
}

template<> inline double _r_factor<false, sixteen>(double b, double invexp)
{
	return sqrt(sqrt(sqrt(sqrt(sqrt(b)))));
}

/*
 * The base component of r_factor
 */
template<bool boxy>
inline double _base(double x, double y, double re, double exponent)
{
	return pow(fabs(x/re), exponent) + pow(fabs(y/re), exponent);
}

template<>
inline double _base<false>(double x, double y, double re, double exponent)
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
inline double _invexp<false>(const double nser, const double exponent)
{
  return nser;
}


/*
 * The main sersic evaluation function for a given X/Y coordinate
 */
template <bool boxy, nser_t t>
inline
double _sersic_for_xy_r(SersicProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	double r_factor;
	if( reuse_r && sp->box == 0. ){
		r_factor = pow(r/sp->re, 1/sp->nser);
	}
	else {
		double base;
		double exponent = sp->box + 2;
		base = _base<boxy>(x, y, sp->re, exponent);
		r_factor = _r_factor<boxy,t>(base,_invexp<boxy>(sp->nser,exponent));
	}

	return exp(-sp->_bn * (r_factor - 1));
}

static inline
void _image_to_sersic_coordinates(SersicProfile *sp, double x, double y, double *x_ser, double *y_ser) {
	x -= sp->xcen;
	y -= sp->ycen;
	*x_ser =  x * sp->_cos_ang + y * sp->_sin_ang;
	*y_ser = -x * sp->_sin_ang + y * sp->_cos_ang;
	*y_ser /= sp->axrat;
}

template <bool boxy, nser_t t>
double _sersic_sumpix(SersicProfile *sp,
                      double x0, double x1, double y0, double y1,
                      unsigned int recur_level, unsigned int max_recursions,
                      unsigned int resolution) {

	double xbin = (x1-x0) / resolution;
	double ybin = (y1-y0) / resolution;
	double half_xbin = xbin/2.;
	double half_ybin = ybin/2.;
	double total = 0, subval, testval;
	double x , y, x_ser, y_ser;
	unsigned int i, j;

	bool recurse = resolution > 1 && recur_level < max_recursions;

	/* The middle X/Y value is used for each pixel */
	x = x0;
	for(i=0; i < resolution; i++) {
		x += half_xbin;
		y = y0;
		for(j=0; j < resolution; j++) {
			y += half_ybin;

			_image_to_sersic_coordinates(sp, x, y, &x_ser, &y_ser);
			subval = _sersic_for_xy_r<boxy, t>(sp, x_ser, y_ser, 0, false);

			if( recurse ) {
				testval = _sersic_for_xy_r<boxy, t>(sp, x_ser, abs(y_ser) + abs(ybin*sp->_cos_ang/sp->axrat), 0, false);
				if( abs(testval/subval - 1.0) > sp->acc ) {
					subval = _sersic_sumpix<boxy, t>(sp,
					                                 x - half_xbin, x + half_xbin,
					                                 y - half_ybin, y + half_ybin,
					                                 recur_level + 1, max_recursions,
					                                 resolution);
				}
			}

			total += subval;
			y += half_ybin;
		}

		x += half_xbin;
	}

	/* Average and return */
	return total / (resolution * resolution);
}

static inline
double sersic_fluxfrac(SersicProfile *sp, double fraction) {
	double ratio = sp->_qgamma(fraction, 2*sp->nser) / sp->_bn;
	return sp->re * pow(ratio, sp->nser);
}

static inline
void sersic_initial_calculations(SersicProfile *sp, Model *model) {

	double nser = sp->nser;
	double re = sp->re;
	double axrat = sp->axrat;
	double mag = sp->mag;
	double box = sp->box + 2;
	double magzero = model->magzero;
	double bn, angrad;

	/*
	 * Calculate the total luminosity used by the sersic profile, used
	 * later to calculate the exact contribution of each pixel.
	 * We save bn back into the profile because it's needed later.
	 */
	sp->_bn = bn = sp->_qgamma(0.5, 2*nser);
	double Rbox = M_PI * box / (4*sp->_beta(1/box, 1 + 1/box));
	double gamma = sp->_gammafn(2*nser);
	double lumtot = pow(re, 2) * 2 * M_PI * nser * gamma * axrat/Rbox * exp(bn)/pow(bn, 2*nser);
	sp->_ie = pow(10, -0.4*(mag - magzero))/lumtot;

	/*
	 * Optionally adjust the user-given re_switch (totally) and resolution
	 * (partially) parameters to more sensible values that will result in faster
	 * profile calculations.
	 */
	if( sp->adjust ) {

		double re_switch;
		unsigned int resolution;

		/*
		 * Find the point at which we capture most of the flux (sensible place
		 * for upscaling). We make sure upscaling doesn't go beyond 20 pixels,
		 * but don't let it become less than 1 pixel (means we do no worse than
		 * GALFIT anywhere)
		 */
		re_switch = ceil(sersic_fluxfrac(sp, 1. - nser*nser/2e3));
		re_switch = max(min(re_switch, 20.), 2.);

		/*
		 * Calculate a bound, adaptive upscale; if re is large then we don't need
		 * so much upscaling
		 */
		resolution = (unsigned int)ceil(160 / re_switch);
		resolution += resolution%2;
		resolution = resolution > 16 ? 16 : resolution;
		resolution = resolution <  4 ?  4 : resolution;

		sp->re_switch = re_switch / re;
		sp->resolution = resolution;

		/*
		 * If the user didn't give a re_max we calculate one that covers
		 * %99.99 of the flux
		 */
		if( sp->re_max == 0 ) {
			sp->re_max = ceil(sersic_fluxfrac(sp, 0.9999));
		}
		sp->_rescale_factor = 1;
		if( sp->rescale_flux ) {
			double flux_r;
			flux_r = bn * pow(sp->re_max/re, 1/nser);
			flux_r = sp->_pgamma(flux_r, 2*nser);
			sp->_rescale_factor = 1/flux_r;
		}

		/* Adjust the accuracy we'll use for sub-pixel integration */
		double acc = 0.4 / nser;
		acc = max(0.1, acc) / axrat;
		sp->acc = acc;

	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into sersic coordinates.
	 *
	 * In galfit the angle started from the Y image axis.
	 */
	angrad = fmod(sp->ang + 90, 360.) * M_PI / 180.;
	sp->_cos_ang = cos(angrad);
	sp->_sin_ang = sin(angrad);

	/* Other way to get sin is doing: cos^2 + sin^2 = 1
	 * sp->_sin_ang = sqrt(1. - cos_ang * cos_ang) * (angrad < M_PI ? -1. : 1.);
	 * The performance seems pretty similar (measured on a x64 Linux with gcc and clang)
	 * and doing sin() is more readable.
	 */

}

/**
 * The sersic validation function
 */
void SersicProfile::validate() {

	if( !this->_pgamma ) {
		throw invalid_parameter("Missing pgamma function on sersic profile");
	}
	if( !this->_qgamma ) {
		throw invalid_parameter("Missing qgamma function on sersic profile");
	}
	if( !this->_gammafn ) {
		throw invalid_parameter("Missing gamma function on sersic profile");
	}
	if( !this->_beta ) {
		throw invalid_parameter("Missing beta function on sersic profile");
	}

}

/**
 * The main sersic evaluation function
 */
template <bool boxy, nser_t t>
void _evaluate(SersicProfile *sp, Model *model, double *image) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_ser, y_ser, r_ser;
	double half_xbin = model->scale_x/2.;
	double half_ybin = model->scale_x/2.;
	double pixel_area = model->scale_x * model->scale_y;

	/*
	 * All the pre-calculations needed by the sersic profile (Ie, cos/sin ang, etc)
	 * We store these profile-global results in the profile structure itself
	 * (it contains extra members to store these values) to avoid passing a long
	 * list of values around every method call.
	 */
	sersic_initial_calculations(sp, model);

	double scale = pixel_area * sp->_ie;
	if( sp->rescale_flux ) {
		scale *= sp->_rescale_factor;
	}

	/* The middle X/Y value is used for each pixel */
	y = 0;
	for(j=0; j < model->height; j++) {
		y += half_ybin;
		x = 0;
		for(i=0; i < model->width; i++) {
			x += half_xbin;

			/* We were instructed to ignore this pixel */
			if( model->calcmask && !model->calcmask[i + j*model->width] ) {
				x += half_xbin;
				continue;
			}

			_image_to_sersic_coordinates(sp, x, y, &x_ser, &y_ser);

			/*
			 * No need for further refinement, return sersic profile
			 * TODO: the radius calculation doesn't take into account boxing
			 */
			r_ser = sqrt(x_ser*x_ser + y_ser*y_ser);
			if( sp->re_max > 0 && r_ser/sp->re > sp->re_max ) {
				pixel_val = 0.;
			}
			else if( sp->rough || r_ser/sp->re > sp->re_switch ){
				pixel_val = _sersic_for_xy_r<boxy, t>(sp, x_ser, y_ser, r_ser, true);
			}
			else {

				bool center = abs(x - sp->xcen) < 1. && abs(y - sp->ycen) < 1.;
				unsigned int resolution = center ? 8 : sp->resolution;
				unsigned int max_recursions = center ? 10 : sp->max_recursions;

				/* Subsample and integrate */
				pixel_val =  _sersic_sumpix<boxy, t>(sp,
				                                     x - half_xbin, x + half_xbin,
				                                     y - half_ybin, y + half_ybin,
				                                     0, max_recursions, resolution);
			}

			image[i + j*model->width] = scale * pixel_val;
			x += half_xbin;
		}
		y += half_ybin;
	}

}

void SersicProfile::evaluate(double *image) {

	Model *m = this->model;

	if( this->box != 0 ) {
		     if( this->nser == 0.5 ) _evaluate<true, pointfive>(this, m, image);
		else if( this->nser == 1 )   _evaluate<true, one>(this, m, image);
		else if( this->nser == 2 )   _evaluate<true, two>(this, m, image);
		else if( this->nser == 3 )   _evaluate<true, three>(this, m, image);
		else if( this->nser == 4 )   _evaluate<true, four>(this, m, image);
		else if( this->nser == 8 )   _evaluate<true, eight>(this, m, image);
		else if( this->nser == 16 )  _evaluate<true, sixteen>(this, m, image);
		else                         _evaluate<true, general>(this, m, image);
	}
	else {
		     if( this->nser == 0.5 ) _evaluate<false, pointfive>(this, m, image);
		else if( this->nser == 1 )   _evaluate<false, one>(this, m, image);
		else if( this->nser == 2 )   _evaluate<false, two>(this, m, image);
		else if( this->nser == 3 )   _evaluate<false, three>(this, m, image);
		else if( this->nser == 4 )   _evaluate<false, four>(this, m, image);
		else if( this->nser == 8 )   _evaluate<false, eight>(this, m, image);
		else if( this->nser == 16 )  _evaluate<false, sixteen>(this, m, image);
		else                         _evaluate<false, general>(this, m, image);
	}
}

#if defined(HAVE_GSL)
double _gsl_qgamma_wrapper(double p, double shape) {
	return gsl_cdf_gamma_Pinv(p, shape, 1);
}
double _gsl_pgamma_wrapper(double q, double shape) {
	return gsl_cdf_gamma_P(q, shape, 1);
}
#elif defined(HAVE_R)
double _Rf_qgamma_wrapper(double p, double shape) {
	return Rf_qgamma(p, shape, 1, 1, 0);
}
double _Rf_pgamma_wrapper(double q, double shape) {
	return Rf_pgamma(q, shape, 1, 1, 0);
}
#endif

/**
 * The sersic creation function
 */
SersicProfile::SersicProfile() :
	Profile()
{

	SersicProfile *p = this;

	/* Sane defaults */
	p->xcen = 0;
	p->ycen = 0;
	p->mag = 15;
	p->re = 1;
	p->nser = 1;
	p->box = 0;
	p->ang   = 0.0;
	p->axrat = 1.;
	p->rough = false;

	p->acc = 0.1;
	p->re_switch = 1.;
	p->resolution = 9;
	p->max_recursions = 2;
	p->adjust = true;

	p->re_max = 0;
	p->rescale_flux = false;

	/*
	 * Point to the corresponding implementation, or leave as NULL if not
	 * possible. In that case the user will have to provide their own functions.
	 */
#if defined(HAVE_GSL)
	p->_qgamma  = &_gsl_qgamma_wrapper;
	p->_pgamma  = &_gsl_pgamma_wrapper;
	p->_gammafn = &gsl_sf_gamma;
	p->_beta    = &gsl_sf_beta;
#elif defined(HAVE_R)
	p->_qgamma  = &_Rf_qgamma_wrapper;
	p->_pgamma  = &_Rf_pgamma_wrapper;
	p->_gammafn = &Rf_gammafn;
	p->_beta    = &Rf_beta;
#else
	p->_qgamma = NULL;
	p->_pgamma = NULL;
	p->_gammafn = NULL;
	p->_beta = NULL;
#endif

}

} /* namespace profit */
