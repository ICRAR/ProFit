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
#include "profit/null.h"
#include "profit/psf.h"
#include "profit/sersic.h"
#include "profit/sky.h"
#include "profit/utils.h"


namespace profit {

Point Model::NO_OFFSET;

Model::Model(unsigned int width, unsigned int height) :
	requested_dimensions(width, height),
	finesampling(1),
	scale(1, 1),
	magzero(0),
	psf(),
	psf_scale(1, 1),
	mask(),
	convolver(),
	crop(true),
	dry_run(false),
	return_finesampled(true),
	opencl_env(),
	omp_threads(0),
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
	if (profile_name == "null") {
		profile = make_shared<NullProfile>(*this, profile_name);
	}
	else if( profile_name == "sky" ) {
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

static
void inform_offset(const Point &offset, Point &offset_out) {
	if (&offset_out != &Model::NO_OFFSET) {
		offset_out = offset;
	}
}

Image Model::evaluate(Point &offset_out) {

	/* Check limits */
	if( !requested_dimensions ) {
		throw invalid_parameter( "Model's requested dimensions are 0");
	}
	else if( this->scale.first <= 0 ) {
		throw invalid_parameter("Model's scale_x cannot be negative or zero");
	}
	else if( this->scale.second <= 0 ) {
		throw invalid_parameter("Model's scale_y cannot be negative or zero");
	}

	/*
	 * If at least one profile is requesting convolving we require
	 * a valid psf.
	 */
	for(auto &profile: this->profiles) {
		if( profile->do_convolve() ) {
			if( !this->psf ) {
				std::ostringstream ss;
				ss << "Profile " << profile->get_name() << " requires convolution but no valid psf was provided";
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
	auto image_dims = requested_dimensions * finesampling;
	Image image(image_dims);
	Point offset(0, 0);

	/* so long folks! */
	if( dry_run ) {
		inform_offset(offset, offset_out);
		return image;
	}

	/*
	 * Generate a separate image for each profile.
	 */
	std::vector<Image> profile_images;
	for(auto &profile: this->profiles) {
		Image profile_image(image_dims);
		profile->adjust_for_finesampling(finesampling);
		profile->evaluate(profile_image, mask, {scale.first / finesampling, scale.second / finesampling}, magzero);
		profile_images.push_back(std::move(profile_image));
	}

	/*
	 * Sum up all results
	 */

	// We first sum up all images that need convolving
	bool do_convolve = false;
	auto it = profile_images.begin();
	for(auto &profile: this->profiles) {
		if( profile->do_convolve() ) {
			do_convolve = true;
			image += *it;
		}
		it++;
	}

	// Now perform convolution on these images
	// The convolution process might produce a larger image;
	// thus we keep track of this bigger size, and the offset
	// of the original image with respect to the new, larger one
	Dimensions conv_dims = image.getDimensions();
	if( do_convolve ) {
		Image psf_img(psf);
		psf_img.normalize();
		if (!convolver) {
			convolver = create_convolver(BRUTE);
		}
		image = convolver->convolve(image, psf_img, mask, crop, offset);
		conv_dims = image.getDimensions();
	}

	// Sum images of profiles that do not require convolution
	Image no_convolved_images(image_dims);
	it = profile_images.begin();
	for(auto &profile: this->profiles) {
		if( !profile->do_convolve() ) {
			no_convolved_images += *it;
		}
		it++;
	}

	// Add non-convolved images on top of the convolved ones
	// taking into account the convolution extension/offset, if any
	if (conv_dims != image_dims) {
		image += no_convolved_images.extend(conv_dims, offset);
	}
	else {
		image += no_convolved_images;
	}

	// Downsample image if necessary
	if (finesampling > 1 && !return_finesampled) {
		image = image.downsample(finesampling, Image::DownsamplingMode::SUM);
		offset /= finesampling;
	}

	/* Done! Good job :-) */
	inform_offset(offset, offset_out);
	return image;
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