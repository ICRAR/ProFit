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

/*
 * We use either GSL or Rmath to provide the low-level
 * beta, gamma and qgamma_inv functions needed by the ferrer profile.
 * If neither is given, the user will have to feed the profiles with
 * the appropriate function pointers after creating them.
 */
#if defined(HAVE_GSL)
#include <gsl/gsl_sf_gamma.h>
#elif defined(HAVE_R)
#define R_NO_REMAP
#include <Rmath.h>
#endif

#include "ferrer.h"

using namespace std;

namespace profit
{

/*
 * The evaluation of the ferrer profile at ferrer coordinates (x,y).
 *
 * The ferrer profile has this form:
 *
 * (1-r_factor)^(a)
 *
 * where r_factor = (r/re)^(2-b)
 *              r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *              B = box parameter
 */


/*
 * The main ferrer evaluation function for a given X/Y coordinate
 */
static inline
double _ferrer_for_xy_r(FerrerProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	double r_factor;
	if( reuse_r && sp->box == 0 ) {
		r_factor = r;
	}
	else if( sp->box == 0 ) {
		r_factor = sqrt(x*x + y*y);
	}
	else {
		double box = sp->box + 2.;
		r_factor = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}

	r_factor /= sp->rout;
	if( r_factor < 1 ) {
		return pow(1 - pow(r_factor, 2 - sp->b), sp->a);
	}
	else {
		return 0;
	}
}

static inline
void _image_to_ferrer_coordinates(FerrerProfile *sp, double x, double y, double *x_ser, double *y_ser) {
	x -= sp->xcen;
	y -= sp->ycen;
	*x_ser =  x * sp->_cos_ang + y * sp->_sin_ang;
	*y_ser = -x * sp->_sin_ang + y * sp->_cos_ang;
	*y_ser /= sp->axrat;
}

double _ferrer_sumpix(FerrerProfile *sp,
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

			_image_to_ferrer_coordinates(sp, x, y, &x_ser, &y_ser);
			subval = _ferrer_for_xy_r(sp, x_ser, y_ser, 0, false);

			if( recurse ) {
				double delta_y_ser = (-xbin*sp->_sin_ang + ybin*sp->_cos_ang)/sp->axrat;
				testval = _ferrer_for_xy_r(sp, abs(x_ser), abs(y_ser) + abs(delta_y_ser), 0, false);
				if( abs(testval/subval - 1.0) > sp->acc ) {
					subval = _ferrer_sumpix(sp,
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
void ferrer_initial_calculations(FerrerProfile *sp, Model *model) {

	double rout = sp->rout;
	double a = sp->a;
	double b = sp->b;
	double axrat = sp->axrat;
	double mag = sp->mag;
	double box = sp->box + 2;
	double magzero = model->magzero;
	double angrad;

	/*
	 * Calculate the total luminosity used by the ferrer profile, used
	 * later to calculate the exact contribution of each pixel.
	 * We save bn back into the profile because it's needed later.
	 */

	/*
	 * lumtot comes from Wolfram Alpha: "integrate between 0 and 1 x(1-x^c)^(a)"
	 * replacement was then made using c=2-b
	 * Further scaling by 2*pi*rout^2 required.
	 * Further scaling by axrat and Rbox.
	 */
	double Rbox = M_PI * box / (4*sp->_beta(1/box, 1 + 1/box));
	double lumtot = pow(rout, 2) * M_PI * (sp->_gammafn(a+1)*sp->_gammafn((4-b)/(2-b))/
	                (sp->_gammafn(sp->a+2/(2-b)+1))) * axrat/Rbox;
	sp->_ie = pow(10, -0.4*(mag - magzero))/lumtot;

	/*
	 * Optionally adjust the user-given re_switch (totally) and resolution
	 * (partially) parameters to more sensible values that will result in faster
	 * profile calculations.
	 */
	if( sp->adjust ) {

		unsigned int resolution;
		/*
		 * Calculate a bound, adaptive upscale; if re is large then we don't need
		 * so much upscaling
		 */
		resolution = (unsigned int)ceil(160 / sp->rout);
		resolution += resolution%2;
		resolution = resolution > 16 ? 16 : resolution;
		resolution = resolution <  4 ?  4 : resolution;

		sp->re_switch = sp->rout;
		sp->resolution = resolution;

		/*
		 * If the user didn't give a re_max we calculate one that covers
		 * %99.99 of the flux
		 */
		if( sp->re_max == 0 ) {
			sp->re_max = sp->rout;
		}

		/* Adjust the accuracy we'll use for sub-pixel integration */
		sp->acc = 0.1/axrat;

	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into ferrer coordinates.
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
 * The ferrer validation function
 */
void FerrerProfile::validate() {

if( !this->_beta ) {
		throw invalid_parameter("Missing beta function on ferrer profile");
	}

}

/**
 * The main ferrer evaluation function
 */
void FerrerProfile::evaluate(double *image) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_ser, y_ser, r_ser;
	double half_xbin = model->scale_x/2.;
	double half_ybin = model->scale_x/2.;
	double pixel_area = model->scale_x * model->scale_y;

	FerrerProfile *sp = this;

	/*
	 * All the pre-calculations needed by the ferrer profile (Ie, cos/sin ang, etc)
	 * We store these profile-global results in the profile structure itself
	 * (it contains extra members to store these values) to avoid passing a long
	 * list of values around every method call.
	 */
	ferrer_initial_calculations(sp, model);

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

			_image_to_ferrer_coordinates(sp, x, y, &x_ser, &y_ser);

			/*
			 * No need for further refinement, return ferrer profile
			 * TODO: the radius calculation doesn't take into account boxing
			 */
			r_ser = sqrt(x_ser*x_ser + y_ser*y_ser);
			if(r_ser > sp->re_max ) {
				pixel_val = 0.;
			}
			else if( sp->rough || r_ser/sp->rout > sp->re_switch ){
				pixel_val = _ferrer_for_xy_r(sp, x_ser, y_ser, r_ser, true);
			}
			else {

				bool center = abs(x - sp->xcen) < 1. && abs(y - sp->ycen) < 1.;
				unsigned int resolution = center ? 8 : sp->resolution;
				unsigned int max_recursions = center ? 10 : sp->max_recursions;

				/* Subsample and integrate */
				pixel_val =  _ferrer_sumpix(sp,
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
 * The ferrer creation function
 */
FerrerProfile::FerrerProfile() :
	Profile()
{

	FerrerProfile *p = this;

	/* Sane defaults */
	p->xcen = 0;
	p->ycen = 0;
	p->mag = 15;
	p->rout= 3;
	p->a= 1;
	p->b= 1;
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
	p->_gammafn = &gsl_sf_gamma;
	p->_beta    = &gsl_sf_beta;
#elif defined(HAVE_R)
	p->_gammafn = &Rf_gammafn;
	p->_beta    = &Rf_beta;
#else
	p->_gammafn = NULL;
	p->_beta = NULL;
#endif

}

} /* namespace profit */
