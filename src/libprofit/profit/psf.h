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

#include "profit/profit.h"

namespace profit
{

/**
 * A PSF profile.
 *
 * PSF profiles simply add the normalized PSF image (for a given magnitude)
 * in a given position onto the model's image.
 */
class PsfProfile : public Profile {

public:

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 */
	PsfProfile(const Model &);

	/*
	 * ---------------------------------------------
	 * Pure virtual functions implementations follow
	 * ---------------------------------------------
	 */
	void validate() override;
	void evaluate(std::vector<double> &image) override;

	/*
	 * -------------------------
	 * Profile parameters follow
	 * -------------------------
	 */

	/**
	 * The X center of this profile
	 */
	double xcen;

	/**
	 * The Y center of this profile
	 */
	double ycen;

	/**
	 * The magnitude of this profile, based on the model's magnitude
	 */
	double mag;

};

} /* namespace profit */

#endif /* _PSF_H_ */
