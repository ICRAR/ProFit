/**
 * Model class implementation
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

#include <sstream>

#include "profit/common.h"
#include "profit/brokenexponential.h"
#include "profit/convolve.h"
#include "profit/coresersic.h"
#include "profit/exceptions.h"
#include "profit/ferrer.h"
#include "profit/king.h"
#include "profit/model.h"
#include "profit/moffat.h"
#include "profit/psf.h"
#include "profit/sersic.h"
#include "profit/sky.h"
#include "profit/utils.h"


namespace profit {

Model::Model(unsigned int width, unsigned int height) :
	width(width), height(height),
	scale_x(1), scale_y(1),
	magzero(0),
	psf(), psf_width(0), psf_height(0),
	psf_scale_x(1), psf_scale_y(1),
	calcmask(),
	convolver(),
	dry_run(false),
#ifdef PROFIT_OPENCL
	opencl_env(),
#endif /* PROFIT_OPENCL */
#ifdef PROFIT_OPENMP
	omp_threads(0),
#endif /* PROFIT_OPENMP */
	profiles()
{
	// no-op
}

bool Model::has_profiles() const {
	return this->profiles.size() > 0;
}

std::shared_ptr<Profile> Model::add_profile(const std::string &profile_name) {

	using std::make_shared;

	std::shared_ptr<Profile> profile;
	if( profile_name == "sky" ) {
		profile = make_shared<SkyProfile>(*this, profile_name);
	}
	else if ( profile_name == "sersic" ) {
		profile = make_shared<SersicProfile>(*this, profile_name);
	}
	else if ( profile_name == "moffat" ) {
		profile = make_shared<MoffatProfile>(*this, profile_name);
	}
	else if ( profile_name == "ferrer" || profile_name == "ferrers" ) {
		profile = make_shared<FerrerProfile>(*this, profile_name);
	}
	else if ( profile_name == "coresersic" ) {
		profile = make_shared<CoreSersicProfile>(*this, profile_name);
	}
	else if ( profile_name == "king" ) {
		profile = make_shared<KingProfile>(*this, profile_name);
	}
	else if ( profile_name == "brokenexp" ) {
		profile = make_shared<BrokenExponentialProfile>(*this, profile_name);
	}
	else if ( profile_name == "psf" ) {
		profile = make_shared<PsfProfile>(*this, profile_name);
	}
	else {
		std::ostringstream ss;
		ss << "Unknown profile name: " << profile_name;
		throw invalid_parameter(ss.str());
	}

	this->profiles.push_back(profile);
	return profile;
}

std::vector<double> Model::evaluate() {

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
		if( profile->do_convolve() ) {
			if( this->psf.empty() ) {
				std::ostringstream ss;
				ss << "Profile " << profile->get_name() << " requires convolution but no psf was provided";
				throw invalid_parameter(ss.str());
			}
			if( !this->psf_width ) {
				throw invalid_parameter("Model's psf width is 0");
			}
			if( !this->psf_height ) {
				throw invalid_parameter("Model's psf height is 0");
			}
			if( this->psf_width * this->psf_height != this->psf.size() ) {
				std::ostringstream ss;
				ss << "PSF dimensions (" << psf_width << "x" << psf_height <<
				      ") don't correspond to PSF length (" << psf.size() << ")";
				throw invalid_parameter(ss.str());
			}
			break;
		}
	}

	/*
	 * Validate all profiles.
	 * Each profile can fail during validation in which case we don't proceed any further
	 */
	for(auto profile: this->profiles) {
		profile->validate();
	}

	// The image we'll eventually return
	Image image(width, height);

	/* so long folks! */
	if( dry_run ) {
		return image.getData();
	}

	/*
	 * Generate a separate image for each profile.
	 */
	std::vector<Image> profile_images;
	for(auto &profile: this->profiles) {
		Image profile_image(width, height);
		profile->evaluate(profile_image.getData());
		profile_images.push_back(std::move(profile_image));
	}

	/*
	 * Sum up all results
	 *
	 * We first sum up all images that need convolving, we convolve them
	 * and after that we add up the remaining images.
	 */
	bool do_convolve = false;
	auto it = profile_images.begin();
	for(auto &profile: this->profiles) {
		if( profile->do_convolve() ) {
			do_convolve = true;
			image += *it;
		}
		it++;
	}
	if( do_convolve ) {
		Image psf_img(psf, psf_width, psf_height);
		psf_img.normalize();
		if (!convolver) {
			convolver = create_convolver(BRUTE);
		}
		Mask mask_img;
		if (!calcmask.empty()) {
			mask_img = Mask(calcmask, width, height);
		}
		image = convolver->convolve(image, psf_img, mask_img);
	}
	it = profile_images.begin();
	for(auto &profile: this->profiles) {
		if( !profile->do_convolve() ) {
			image += *it;
		}
		it++;
	}

	/* Done! Good job :-) */
	return image.getData();
}

std::map<std::string, std::shared_ptr<ProfileStats>> Model::get_stats() const {
	std::map<std::string, std::shared_ptr<ProfileStats>> stats;
	for(auto p: profiles) {
		stats[p->get_name()] = p->get_stats();
	}
	return stats;
}

#ifdef PROFIT_DEBUG
std::map<std::string, std::map<int, int>> Model::get_profile_integrations() const {
	std::map<std::string, std::map<int, int>> profile_integrations;
	for(auto p: profiles) {
		RadialProfile *rp = dynamic_cast<RadialProfile *>(p.get());
		if( rp ) {
			profile_integrations[rp->get_name()] = rp->get_integrations();
		}
	}
	return profile_integrations;
}
#endif /* PROFIT_DEBUG */

} /* namespace profit */