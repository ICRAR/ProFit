/**
 * Header file for PSF profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
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
#ifndef _PSF_H_
#define _PSF_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include "profit.h"

typedef struct _profit_psf_profile {
	profit_profile profile;
	double xcen;
	double ycen;
	double mag;

	/* This is calculated from mag */
	double scale;
} profit_psf_profile;

profit_profile *profit_create_psf(void);

#ifdef __cplusplus
}
#endif

#endif /* _PSF_H_ */
