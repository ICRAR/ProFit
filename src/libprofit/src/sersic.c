/**
 * Sersic profile implementation
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

#include <math.h>
#include <stdlib.h>
#include <string.h>

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

static inline
double _sersic_for_xy_r(profit_sersic_profile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

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
	 *  r_factor = (x/re)^{2+b} + (y/re)^{2+b})^{1/(nser*(2+b)}
	 *
	 */

	double r_factor;
	if( reuse_r && sp->box == 0. ){
		r_factor = pow(r/sp->re, 1/sp->nser);
	}
	else {
		double exponent = sp->box + 2;
		double base = pow(fabs(x/sp->re), exponent) + pow(fabs(y/sp->re), exponent);
		double divisor = sp->nser*exponent;

		if( divisor == 0.5 ) {
			r_factor = base*base;
		} else if( divisor == 1. ) {
			r_factor = base;
		} else if( divisor == 2. ) {
			r_factor = sqrt(base);
		} else if( divisor == 3. ) {
			r_factor = cbrt(base);
		} else if( divisor == 4. ) {
			r_factor = sqrt(sqrt((base)));
		} else {
			r_factor = pow(base, 1/divisor);
		}
	}

	return exp(-sp->_bn * (r_factor - 1));
}

static inline
void _image_to_sersic_coordinates(profit_sersic_profile *sp, double x, double y, double *x_ser, double *y_ser) {
	x -= sp->xcen;
	y -= sp->ycen;
	*x_ser =  x * sp->_cos_ang + y * sp->_sin_ang;
	*y_ser = -x * sp->_sin_ang + y * sp->_cos_ang;
	*y_ser /= sp->axrat;
}

static
double _sersic_sumpix(profit_sersic_profile *sp,
                      double x0, double x1, double y0, double y1,
                      unsigned int recur_level) {

	double xbin = (x1-x0) / sp->resolution;
	double ybin = (y1-y0) / sp->resolution;
	double half_xbin = xbin/2.;
	double half_ybin = ybin/2.;
	double total = 0, subval, testval;
	double x , y, x_ser, y_ser;
	unsigned int i, j;

	bool recurse = sp->resolution > 1 && recur_level < sp->max_recursions;

	/* The middle X/Y value is used for each pixel */
	x = x0;
	for(i=0; i < sp->resolution; i++) {
		x += half_xbin;
		y = y0;
		for(j=0; j < sp->resolution; j++) {
			y += half_ybin;

			_image_to_sersic_coordinates(sp, x, y, &x_ser, &y_ser);
			subval = _sersic_for_xy_r(sp, x_ser, y_ser, 0, false);

			if( recurse ) {
				testval = _sersic_for_xy_r(sp, x_ser, fabs(y_ser) + fabs(ybin/sp->axrat), 0, false);
				if( fabs(testval/subval - 1.0) > sp->acc ) {
					subval = _sersic_sumpix(sp,
					                        x - half_xbin, x + half_xbin,
					                        y - half_ybin, y + half_ybin,
					                        recur_level + 1);
				}
			}

			total += subval;
			y += half_ybin;
		}

		x += half_xbin;
	}

	/* Average and return */
	return total / (sp->resolution * sp->resolution);
}

static
void profit_make_sersic(profit_profile *profile, profit_model *model, double *image) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_ser, y_ser, r_ser;
	double half_xbin = model->xbin/2.;
	double half_ybin = model->ybin/2.;
	double bin_area = model->xbin * model->ybin;
	profit_sersic_profile *sp = (profit_sersic_profile *)profile;

	/* The middle X/Y value is used for each pixel */
	x = 0;
	for(i=0; i < model->width; i++) {
		x += half_xbin;
		y = 0;
		for(j=0; j < model->height; j++) {
			y += half_ybin;

			_image_to_sersic_coordinates(sp, x, y, &x_ser, &y_ser);

			/*
			 * No need for further refinement, return sersic profile
			 * TODO: the radius calculation doesn't take into account boxing
			 */
			r_ser = sqrt(x_ser*x_ser + y_ser*y_ser);
			if( sp->rough || r_ser/sp->re > sp->re_switch ){
				pixel_val = _sersic_for_xy_r(sp, x_ser, y_ser, r_ser, true);
			}
			else {
				/* Subsample and integrate */
				pixel_val =  _sersic_sumpix(sp,
				                            x - model->xbin/2, x + model->xbin/2,
				                            y - model->ybin/2, y + model->ybin/2,
				                            0);
			}

			image[i + j*model->width] = bin_area * sp->_ie * pixel_val;
			y += half_ybin;
		}
		x += half_xbin;
	}

}

#include <stdio.h>
void dump_profile(profit_sersic_profile *s) {
	printf("Sersic profile details:\n");
	printf("Shape paremeters:\n");
	printf("  xcen = %f\n", s->xcen);
	printf("  ycen = %f\n", s->ycen);
	printf("   mag = %f\n", s->mag);
	printf("    re = %f\n", s->re);
	printf("  nser = %f\n", s->nser);
	printf("   box = %f\n", s->box);
	printf("   ang = %f\n", s->ang);
	printf(" axrat = %f\n", s->axrat);
	printf("Sub-pixel integration parameters:\n");
	printf("          rough = %d\n", s->rough);
	printf("            acc = %f\n", s->xcen);
	printf("      re_switch = %f\n", s->re_switch);
	printf("     resolution = %u\n", s->resolution);
	printf(" max_recursions = %u\n", s->max_recursions);
}

static
void profit_init_sersic(profit_profile *profile, profit_model *model) {

	profit_sersic_profile *sersic_p = (profit_sersic_profile *)profile;
	double nser = sersic_p->nser;
	double re = sersic_p->re;
	double axrat = sersic_p->axrat;
	double mag = sersic_p->mag;
	double box = sersic_p->box + 2;
	double magzero = model->magzero;
	double bn, angrad;

	if( !sersic_p->_qgamma ) {
		profile->error = strdup("Missing qgamma function on sersic profile");
		return;
	}
	if( !sersic_p->_gammafn ) {
		profile->error = strdup("Missing gamma function on sersic profile");
		return;
	}
	if( !sersic_p->_beta ) {
		profile->error = strdup("Missing beta function on sersic profile");
		return;
	}

	/*
	 * Calculate the total luminosity used by the sersic profile, used
	 * later to calculate the exact contribution of each pixel.
	 * We save bn back into the profile because it's needed later.
	 */
	sersic_p->_bn = bn = sersic_p->_qgamma(0.5, 2*nser, 1);
	double Rbox = M_PI * box / (4*sersic_p->_beta(1/box, 1 + 1/box));
	double gamma = sersic_p->_gammafn(2*nser);
	double lumtot = pow(re, 2) * 2 * M_PI * nser * gamma * axrat/Rbox * exp(bn)/pow(bn, 2*nser);
	sersic_p->_ie = pow(10, -0.4*(mag - magzero))/lumtot;

	/*
	 * Optionally adjust the user-given re_switch (totally) and resolution
	 * (partially) parameters to more sensible values that will result in faster
	 * profile calculations.
	 */
	if( sersic_p->adjust ) {

		double re_switch, flux_frac;
		unsigned int resolution;

		/*
		 * Find the point at which we capture most of the flux (sensible place
		 * for upscaling). We make sure upscaling doesn't go beyond 20 pixels,
		 * but don't let it become less than 1 pixel (means we do no worse than
		 * GALFIT anywhere)
		 */
		flux_frac = 1 - (nser*nser)/1e3;
		re_switch = ceil( re * pow(sersic_p->_qgamma(flux_frac, 2*nser, 1) / bn, nser) );
		re_switch = fmax(fmin(re_switch, 20.), 1);
		re_switch /= re;

		/*
		 * Calculate a bound, adaptive upscale; if re is large then we don't need
		 * so much upscaling
		 */
		resolution = (unsigned int)ceil(sersic_p->resolution * sersic_p->resolution / re_switch);
		resolution = resolution > 9 ? 9 : resolution;
		resolution = resolution < 3 ? 3 : resolution;

		sersic_p->re_switch = re_switch;
		sersic_p->resolution = resolution;
	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into sersic coordinates.
	 *
	 * In galfit the angle started from the Y image axis.
	 */
	angrad = fmod(sersic_p->ang + 90, 360.) * M_PI / 180.;
	sersic_p->_cos_ang = cos(angrad);
	sersic_p->_sin_ang = sin(angrad);

	/* Other way to get sin is doing: cos^2 + sin^2 = 1
	 * sersic_p->_sin_ang = sqrt(1. - cos_ang * cos_ang) * (angrad < M_PI ? -1. : 1.);
	 * The performance seems pretty similar (measured on a x64 Linux with gcc and clang)
	 * and doing sin() is more readable.
	 */

	//dump_profile(sersic_p);

}

#if defined(HAVE_R)
double _Rf_qgamma_wrapper(double a, double b, double c) {
	return Rf_qgamma(a, b, c, 1, 0);
}
#endif

profit_profile *profit_create_sersic() {
	profit_sersic_profile *p = (profit_sersic_profile *)malloc(sizeof(profit_sersic_profile));
	p->profile.init_profile = &profit_init_sersic;
	p->profile.make_profile = &profit_make_sersic;

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

	/*
	 * Point to the corresponding implementation, or leave as NULL if not
	 * possible. In that case the user will have to provide their own functions.
	 */
#if defined(HAVE_GSL)
	p->_beta = &gsl_sf_beta;
	p->_gammafn = &gsl_sf_gamma;
	p->_qgamma = &gsl_cdf_gamma_Qinv;
#elif defined(HAVE_R)
	p->_qgamma = &_Rf_qgamma_wrapper;
	p->_gammafn = &Rf_gammafn;
	p->_beta = &Rf_beta;
#else
	p->_qgamma = NULL;
	p->_gammafn = NULL;
	p->_beta = NULL;
#endif
	return (profit_profile *)p;
}
