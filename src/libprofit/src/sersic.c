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

#include "sersic.h"

static
double profit_sumpix(double xcen, double ycen, double xlim0, double xlim1, double ylim0, double ylim1,
                     double re, double nser, double angrad, double axrat, double box,
                     double bn, int N, int recur, double acc) {

	double rad,x,y,xmod,ymod,radmod,angmod;
	double xbin=(xlim1-xlim0)/N;
	double ybin=(ylim1-ylim0)/N;
	double sumpixel=0, addval, prevaddval = 0;
	int upscale=20;
	x=xlim0;
	for(int i = 0; i < N; i++) {
		recur=0;
		y=ylim0;
		for(int j = 0; j < N; j++) {
			rad=sqrt(pow(x+xbin/2-xcen,2)+pow(y+ybin/2-ycen,2));
			angmod=atan2(x+xbin/2-xcen,y+ybin/2-ycen)-angrad;
			xmod=rad*sin(angmod);
			ymod=rad*cos(angmod);
			xmod=xmod/axrat;
			radmod=pow(pow(fabs(xmod),box)+pow(fabs(ymod),box),1/(box));
			addval=exp(-bn*(pow(radmod/re,1/nser)-1));
			if( j > 0 && recur < 3 ){
				if( (addval/prevaddval > (1+acc)) || (addval/prevaddval < 1/(1+acc)) ){
					recur++;
					addval=profit_sumpix(xcen,ycen,x, x+xbin, y, y+ybin,re,nser,angrad,axrat,box,bn,upscale,recur,acc);
				}
			}
			sumpixel+=addval;
			prevaddval=addval;
			y=y+ybin;
		}
		x=x+xbin;
	}

	return(sumpixel/pow(N,2));
}

static
int _sersic_at_xy(profit_sersic_profile *sp,
                  profit_model *model,
                  unsigned int x, unsigned int y,
                  double *result) {

	double angrad = -sp->ang * M_PI/180;
	double re   = sp->re;
	double nser = sp->nser;
	double box  = sp->box + 2;
	double xbin = model->xbin;
	double ybin = model->ybin;
	/* unsigned int depth = 0; */

	/*
	 * Transform the X/Y position so it accounts for translation from the
	 * profile center, rotation and ellipse scaling in the X axis
	 */
	double rad = sqrt( pow(x+xbin/2-sp->xcen,2) + pow(y+ybin/2-sp->ycen,2) );
	double angmod = atan2(x-sp->xcen, y-sp->ycen) - angrad;
	double xmod = rad * sin(angmod) / sp->axrat;
	double ymod = rad * cos(angmod);
	double radmod = pow( pow(fabs(xmod),box) + pow(fabs(ymod),box), 1/(box));

	unsigned int upscale = 4;

	/*
	 * No need for further refinement, return sersic profile
	 */
	if( sp->rough || radmod > 2*re ){
		*result = exp( -sp->bn * (pow(radmod/re, 1/nser) - 1) );
		*result *= xbin*ybin*sp->Ie;
		return 0;
	}

	/*
	 * Adaptive scaling: perform subsampling and return sum of subsampling
	 */
	double locscale = xbin/radmod;
	locscale = (locscale > 10 ) ? 10 : locscale;
	if(radmod<xbin) {
		upscale = ceil(8*nser*locscale);
		/*depth=3;*/
	}
	else if( radmod < 0.1*re ){
		upscale = ceil(8*nser*locscale);
		/*depth=2;*/
	}
	else if( radmod < 0.25*re ){
		upscale = ceil(4*nser*locscale);
		/*depth=2;*/
	}
	else if( radmod < 0.5*re ){
		upscale = ceil(2*nser*locscale);
		/*depth=2;*/
	}
	else if( radmod < re ){
		upscale = ceil(nser*locscale);
		/*depth=1;*/
	}
	else if( radmod <= 2*re ){
		upscale = ceil((nser/2)*locscale);
		/*depth=0;*/
	}

	/* Min/max scale */
	upscale = (upscale < 4) ? 4 : upscale;
	upscale = (upscale > 100) ? 100 : upscale;

	*result = profit_sumpix(sp->xcen, sp->ycen, x, x+xbin, y, y+ybin,
	                        re, nser, angrad, sp->axrat, box, sp->bn,
	                        upscale, 0, 1e-1);
	*result *= xbin * ybin * sp->Ie;
	return 0;

}

void profit_make_sersic(profit_profile *profile, profit_model *model, double *image) {

	unsigned int i, j;
	double x, y;
	profit_sersic_profile *sp = (profit_sersic_profile *)profile;

	for(i=0; i < model->width; i++) {
		x = model->xbin * i;
		for(j=0; j < model->height; j++) {
			y = model->ybin * j;
			double pix_val;
			_sersic_at_xy(sp, model, x, y, &pix_val);
			image[j*model->width + i] = pix_val;
		}
	}

}

void profit_init_sersic(profit_profile *profile, profit_model *model) {

	profit_sersic_profile *sersic_p = (profit_sersic_profile *)profile;
	double nser = sersic_p->nser;
	double re = sersic_p->re;
	double axrat = sersic_p->axrat;
	double mag = sersic_p->mag;
	double box = sersic_p->box + 2;
	double magzero = model->magzero;
	double bn;

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
	 * Calculate the total luminosity needed by the sersic profile
	 * and save it back on the sersic profile
	 */
	sersic_p->bn = bn = sersic_p->_qgamma(0.5, 2*nser, 1);
	double Rbox = M_PI * box / (4*sersic_p->_beta(1/box, 1 + 1/box));
	double gamma = sersic_p->_gammafn(2*nser);
	double lumtot = pow(re, 2) * 2 * M_PI * nser * gamma * axrat/Rbox * exp(bn)/pow(bn, 2*nser);
	sersic_p->Ie = pow(10, -0.4*(mag-magzero))/lumtot;

	return;
}

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
	p->rough = 0;
	p->_qgamma = NULL;
	p->_gammafn = NULL;
	p->_beta = NULL;
	return (profit_profile *)p;
}
