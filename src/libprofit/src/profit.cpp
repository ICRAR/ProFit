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

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "profit/brokenexponential.h"
#include "profit/convolve.h"
#include "profit/coresersic.h"
#include "profit/ferrer.h"
#include "profit/king.h"
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

Profile::Profile(const Model &model) :
	model(model),
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
	psf(), psf_width(0), psf_height(0),
	psf_scale_x(1), psf_scale_y(1),
	calcmask(), profiles()
{
	// no-op
}

bool Model::has_profiles() const {
	return this->profiles.size() > 0;
}

Profile &Model::add_profile(const string &profile_name) {

	Profile * profile = nullptr;
	if( profile_name == "sky" ) {
		profile = static_cast<Profile *>(new SkyProfile(*this));
	}
	else if ( profile_name == "sersic" ) {
		profile = static_cast<Profile *>(new SersicProfile(*this));
	}
	else if ( profile_name == "moffat" ) {
		profile = static_cast<Profile *>(new MoffatProfile(*this));
	}
	else if ( profile_name == "ferrer" || profile_name == "ferrers" ) {
		profile = static_cast<Profile *>(new FerrerProfile(*this));
	}
	else if ( profile_name == "coresersic" ) {
		profile = static_cast<Profile *>(new CoreSersicProfile(*this));
	}
	else if ( profile_name == "king" ) {
		profile = static_cast<Profile *>(new KingProfile(*this));
	}
	else if ( profile_name == "brokenexp" ) {
		profile = static_cast<Profile *>(new BrokenExponentialProfile(*this));
	}
	else if ( profile_name == "psf" ) {
		profile = static_cast<Profile *>(new PsfProfile(*this));
	}
	else {
		ostringstream ss;
		ss << "Unknown profile name: " << profile_name;
		throw invalid_parameter(ss.str());
	}

	profile->name = profile_name;
	this->profiles.push_back(profile);
	return *profile;
}

vector<double> Model::evaluate() {

	/* Check limits */
	if( !this->width ) {
		throw invalid_parameter( "Model's width is 0");
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
			if( this->psf.empty() ) {
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

	vector<double> image(this->width * this->height, 0);

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
	vector<vector<double>> profile_images;
	for(auto profile: this->profiles) {
		vector<double> profile_image(this->width * this->height, 0);
		profile->evaluate(profile_image);
		profile_images.push_back(profile_image);
	}

	/*
	 * Sum up all results
	 *
	 * We first sum up all images that need convolving, we convolve them
	 * and after that we add up the remaining images.
	 */
	bool do_convolve = false;
	auto it = profile_images.begin();
	for(auto profile: this->profiles) {
		if( profile->convolve ) {
			do_convolve = true;
			add_images(image, *it);
		}
		it++;
	}
	if( do_convolve ) {
		vector<double> psf(this->psf);
		normalize(psf);
		image = convolve(image, this->width, this->height, psf, this->psf_width, this->psf_height, this->calcmask);
	}
	it = profile_images.begin();
	for(auto profile: this->profiles) {
		if( !profile->convolve ) {
			add_images(image, *it);
		}
		it++;
	}

	/* Done! Good job :-) */
	return image;
}

Model::~Model() {
	for(auto profile: this->profiles) {
		delete profile;
	}
}

} /* namespace profit */
