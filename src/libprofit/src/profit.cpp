/**
 * libprofit main entry-point routines
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

#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "profit/convolve.h"
#include "profit/ferrer.h"
#include "profit/moffat.h"
#include "profit/profit.h"
#include "profit/psf.h"
#include "profit/sersic.h"
#include "profit/sky.h"
#include "profit/utils.h"

using namespace std;

namespace profit {

invalid_parameter::invalid_parameter(const string &what_arg) :
	exception(),
	m_what(what_arg)
{
	// no-op
}

invalid_parameter::invalid_parameter(const invalid_parameter &e) :
   m_what(e.m_what)
{
	// no-op
}

invalid_parameter::~invalid_parameter() throw () {
	// no-op
}

const char *invalid_parameter::what() const throw() {
	return m_what.c_str();
}

Profile::Profile() :
	convolve(false)
{
	// no-op
}

Profile::~Profile()
{
	// no-op
}

Model::Model() :
	width(0), height(0),
	scale_x(1), scale_y(1),
	magzero(0),
	psf(NULL), psf_width(0), psf_height(0),
	psf_scale_x(1), psf_scale_y(1),
	calcmask(NULL), image(NULL),
	profiles()
{
	// no-op
}

Profile* Model::add_profile(string profile_name) {

	Profile *profile = NULL;
	if( profile_name == "sky" ) {
		profile = static_cast<Profile *>(new SkyProfile());
	}
	else if ( profile_name == "sersic" ) {
		profile = static_cast<Profile *>(new SersicProfile());
	}
	else if ( profile_name == "moffat" ) {
		profile = static_cast<Profile *>(new MoffatProfile());
	}
	else if ( profile_name == "ferrer" ) {
		profile = static_cast<Profile *>(new FerrerProfile());
	}
	else if ( profile_name == "psf" ) {
		profile = static_cast<Profile *>(new PsfProfile());
	}

	if( profile == NULL ) {
		return NULL;
	}

	profile->model = this;
	profile->name = profile_name;
	this->profiles.push_back(profile);
	return profile;
}

void Model::evaluate() {

	/* Check limits */
	if( !this->width ) {
		throw invalid_parameter( "Model's width is 0");
		return;
	}
	else if( !this->height ) {
		throw invalid_parameter("Model's height is 0");
	}
	else if( this->scale_x <= 0 ) {
		throw invalid_parameter("Model's scale_x cannot be negative or zero");
	}
	else if( this->scale_y <= 0 ) {
		throw invalid_parameter("Model's scale_y cannot be negative or zero");
	}

	/*
	 * If at least one profile is requesting convolving we require
	 * a valid psf.
	 */
	for(auto profile: this->profiles) {
		if( profile->convolve ) {
			if( !this->psf ) {
				stringstream ss;
				ss <<  "Profile " << profile->name << " requires convolution but no psf was provided";
				throw invalid_parameter(ss.str());
			}
			if( !this->psf_width ) {
				throw invalid_parameter("Model's psf width is 0");
			}
			if( !this->psf_height ) {
				throw invalid_parameter("Model's psf height is 0");
			}
			break;
		}
	}

	this->image = new double[this->width * this->height];
	memset(this->image, 0, sizeof(double) * this->width * this->height);

	/*
	 * Validate all profiles.
	 * Each profile can fail during validation in which case we don't proceed any further
	 */
	for(auto profile: this->profiles) {
		profile->validate();
	}

	/*
	 * Generate a separate image for each profile.
	 *
	 * We optionally use OpenMP to parallelize this per-profile image
	 * generation. Depending on how many there are we might get a speed up, so
	 * probably we should study what is the best way to go here (e.g.,
	 * parallelize only if we have more than 2 or 3 profiles)
	 */
	vector<double *> profile_images;
	for(auto profile: this->profiles) {
		double *profile_image = new double[this->width * this->height];
		memset(profile_image, 0, sizeof(double) * this->width * this->height);
		profile->evaluate(profile_image);
		profile_images.push_back(profile_image);
	}

	/*
	 * Sum up all results
	 *
	 * We first sum up all images that need convolving, we convolve them
	 * and after that we add up the remaining images.
	 */
	bool convolve = false;
	vector<double *>::iterator it = profile_images.begin();
	for(auto profile: this->profiles) {
		if( profile->convolve ) {
			convolve = true;
			add_images(this->image, *it, this->width, this->height);
		}
		it++;
	}
	if( convolve ) {
		size_t psf_size = sizeof(double) * this->psf_width * this->psf_height;
		double* psf = new double[psf_size];
		memcpy(psf, this->psf, psf_size);
		normalize(psf, this->psf_width, this->psf_height);
		profit::convolve(this->image, this->width, this->height, psf, this->psf_width, this->psf_height, this->calcmask, true);
		delete [] psf;
	}
	it = profile_images.begin();
	for(auto profile: this->profiles) {
		if( !profile->convolve ) {
			add_images(this->image, *it, this->width, this->height);
		}
		delete [] *it;
		it++;
	}

	/* Done! Good job :-) */
}

Model::~Model() {

	for(auto profile: this->profiles) {
		delete profile;
	}
	delete [] this->image;
	delete [] this->psf;
	delete [] this->calcmask;
}

} /* namespace profit */
