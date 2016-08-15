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

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

/*
 * We use either GSL or Rmath to provide the low-level
 * beta, gamma and qgamma_inv functions needed by the moffat profile.
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

#include "moffat.h"

using namespace std;

namespace profit
{

/*
 * The evaluation of the moffat profile at moffat coordinates (x,y).
 *
 * The moffat profile has this form:
 *
 * e^{-bn * (r_factor - 1)}
 
 * where r_factor = (r/Re)^{1/con}
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
 * The main moffat evaluation function for a given X/Y coordinate
 */
inline
double _moffat_for_xy_r(MoffatProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

    double r_factor = (sp->box == 0) ?
               sqrt(x*x + y*y)/sp->_re :
               (pow(pow(abs(x),2.+sp->box)+pow(abs(y),2.+sp->box),1./(2.+sp->box)) ) / sp->_re;

	return 1/(pow(1+pow(r_factor,2), sp->con));
}

static inline
void _image_to_moffat_coordinates(MoffatProfile *sp, double x, double y, double *x_ser, double *y_ser) {
	x -= sp->xcen;
	y -= sp->ycen;
	*x_ser =  x * sp->_cos_ang + y * sp->_sin_ang;
	*y_ser = -x * sp->_sin_ang + y * sp->_cos_ang;
	*y_ser /= sp->axrat;
}

double _moffat_sumpix(MoffatProfile *sp,
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

			_image_to_moffat_coordinates(sp, x, y, &x_ser, &y_ser);
			subval = _moffat_for_xy_r(sp, x_ser, y_ser, 0, false);

			if( recurse ) {
				testval = _moffat_for_xy_r(sp, x_ser, abs(y_ser) + abs(ybin*sp->_cos_ang/sp->axrat), 0, false);
				if( abs(testval/subval - 1.0) > sp->acc ) {
					subval = _moffat_sumpix(sp,
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
double moffat_fluxfrac(MoffatProfile *sp, double fraction) {
    return 1;
}

static inline
void moffat_initial_calculations(MoffatProfile *sp, Model *model) {

	double con = sp->con;
	double fwhm = sp->fwhm;
	double axrat = sp->axrat;
	double mag = sp->mag;
	double box = sp->box + 2;
	double magzero = model->magzero;
	double angrad;

	/*
	 * Calculate the total luminosity used by the moffat profile, used
	 * later to calculate the exact contribution of each pixel.
	 * We save bn back into the profile because it's needed later.
	 */
	double Rbox = M_PI * box / (4*sp->_beta(1/box, 1 + 1/box));
    double re = sp->_re = fwhm/(2*sqrt(pow(2,(1/con))-1));
    double lumtot = pow(re, 2) * 2 * M_PI * axrat/(con-1)/Rbox;
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
		re_switch = ceil(moffat_fluxfrac(sp, 1. - con*con/2e3));
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
			sp->re_max = ceil(moffat_fluxfrac(sp, 0.9999));
		}

		/* Adjust the accuracy we'll use for sub-pixel integration */
		double acc = 0.4 / con;
		acc = max(0.1, acc) / axrat;
		sp->acc = acc;

	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into moffat coordinates.
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
 * The moffat validation function
 */
void MoffatProfile::validate() {

if( !this->_beta ) {
		throw invalid_parameter("Missing beta function on moffat profile");
	}

}

/**
 * The main moffat evaluation function
 */
void MoffatProfile::evaluate(double *image) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_ser, y_ser, r_ser;
	double half_xbin = model->scale_x/2.;
	double half_ybin = model->scale_x/2.;
	double pixel_area = model->scale_x * model->scale_y;
	
	MoffatProfile *sp = this;

    /* We never convolve */
    this->convolve = false;

	/*
	 * All the pre-calculations needed by the moffat profile (Ie, cos/sin ang, etc)
	 * We store these profile-global results in the profile structure itself
	 * (it contains extra members to store these values) to avoid passing a long
	 * list of values around every method call.
	 */
	moffat_initial_calculations(sp, model);

	double scale = pixel_area * sp->_ie;

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

			_image_to_moffat_coordinates(sp, x, y, &x_ser, &y_ser);

			/*
			 * No need for further refinement, return moffat profile
			 * TODO: the radius calculation doesn't take into account boxing
			 */
			r_ser = sqrt(x_ser*x_ser + y_ser*y_ser);
			if( sp->re_max > 0 && r_ser/sp->_re > sp->re_max ) {
				pixel_val = 0.;
			}
			else if( sp->rough || r_ser/sp->_re > sp->re_switch ){
				pixel_val = _moffat_for_xy_r(sp, x_ser, y_ser, r_ser, true);
			}
			else {

				bool center = abs(x - sp->xcen) < 1. && abs(y - sp->ycen) < 1.;
				unsigned int resolution = center ? 8 : sp->resolution;
				unsigned int max_recursions = center ? 10 : sp->max_recursions;

				/* Subsample and integrate */
				pixel_val =  _moffat_sumpix(sp,
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

/**
 * The moffat creation function
 */
MoffatProfile::MoffatProfile() :
	Profile()
{

	MoffatProfile *p = this;

	/* Sane defaults */
	p->xcen = 0;
	p->ycen = 0;
	p->mag = 15;
	p->fwhm = 1;
	p->con= 2;
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

	/*
	 * Point to the corresponding implementation, or leave as NULL if not
	 * possible. In that case the user will have to provide their own functions.
	 */
#if defined(HAVE_GSL)
	p->_beta    = &gsl_sf_beta;
#elif defined(HAVE_R)
	p->_beta    = &Rf_beta;
#else
	p->_beta = NULL;
#endif

}

} /* namespace profit */
