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

#include <algorithm>
#include <cassert>
#include <functional>
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
	adjust_mask(true),
	crop(true),
	dry_run(false),
	return_finesampled(true),
	opencl_env(),
	omp_threads(0),
	profiles()
{
	// no-op
}

Model::Model(Dimensions dimensions) :
	requested_dimensions(dimensions),
	finesampling(1),
	scale(1, 1),
	magzero(0),
	psf(),
	psf_scale(1, 1),
	mask(),
	convolver(),
	adjust_mask(true),
	crop(true),
	dry_run(false),
	return_finesampled(true),
	opencl_env(),
	omp_threads(0),
	profiles()
{
}

bool Model::has_profiles() const {
	return this->profiles.size() > 0;
}

template <typename P>
ProfilePtr Model::make_profile(const std::string &profile_name)
{
	auto profile = std::make_shared<P>(*this, profile_name);
	this->profiles.push_back(profile);
	return profile;
}

ProfilePtr Model::add_profile(const std::string &profile_name)
{
	if (profile_name == "null") {
		return make_profile<NullProfile>(profile_name);
	}
	else if (profile_name == "sky") {
		return make_profile<SkyProfile>(profile_name);
	}
	else if (profile_name == "sersic") {
		return make_profile<SersicProfile>(profile_name);
	}
	else if (profile_name == "moffat") {
		return make_profile<MoffatProfile>(profile_name);
	}
	else if (profile_name == "ferrer" || profile_name == "ferrers") {
		return make_profile<FerrerProfile>(profile_name);
	}
	else if (profile_name == "coresersic") {
		return make_profile<CoreSersicProfile>(profile_name);
	}
	else if (profile_name == "king") {
		return make_profile<KingProfile>(profile_name);
	}
	else if (profile_name == "brokenexp") {
		return make_profile<BrokenExponentialProfile>(profile_name);
	}
	else if (profile_name == "psf") {
		return make_profile<PsfProfile>(profile_name);
	}
	std::ostringstream ss;
	ss << "Unknown profile name: " << profile_name;
	throw invalid_parameter(ss.str());
}

static
void inform_offset(const Point &offset, Point &offset_out) {
	if (&offset_out != &Model::NO_OFFSET) {
		offset_out = offset;
	}
}

ConvolverPtr &Model::ensure_convolver()
{
	if (!convolver) {
		convolver = create_convolver(BRUTE);
	}
	return convolver;
}

void Model::analyze_expansion_requirements(const Dimensions &dimensions,
    const Mask &mask, const Image &psf, unsigned int finesampling,
    input_analysis& analysis, bool adjust_mask)
{
	// In order to correctly capture all flux in the model when convolving
	// we need to produce images that span outside the originally requested
	// area. This way, when convolution happens the pixels in the original
	// region will contain flux coming from "outside" the image.
	// If the mask is conveniently smaller and centrally located over the image
	// then we can actually avoid having to generated bigger model images
	analysis.mask_needs_convolution = false;
	analysis.mask_needs_psf_padding = false;
	bool model_needs_psf_padding = analysis.convolution_required;
	if (mask && !adjust_mask) {
		model_needs_psf_padding = mask.getDimensions() > dimensions * finesampling;
	}
	else if (mask && analysis.convolution_required) {
		auto bounds = mask.bounding_box() * finesampling;
		auto mask_pad_low = bounds.first;
		auto mask_pad_up = mask.getDimensions() * finesampling - bounds.second;
		auto needed = psf.getDimensions() / 2;
		model_needs_psf_padding = !(mask_pad_low >= needed) || !(mask_pad_up >= needed);
		analysis.mask_needs_convolution = true;
		analysis.mask_needs_psf_padding = model_needs_psf_padding;
	}
	if (model_needs_psf_padding) {
		analysis.psf_padding = psf.getDimensions() / 2;
	} else {
		analysis.psf_padding = Dimensions{0, 0};
	}
	analysis.mask_needs_adjustment = needs_adjustment(mask, finesampling, analysis);
	analysis.drawing_dims = dimensions * finesampling + analysis.psf_padding * 2;
}

Model::input_analysis Model::analyze_inputs() const
{
	/* Check limits */
	if (!requested_dimensions) {
		throw invalid_parameter( "Model's requested dimensions are 0");
	}
	else if (scale.first <= 0) {
		throw invalid_parameter("Model's scale_x cannot be negative or zero");
	}
	else if (scale.second <= 0) {
		throw invalid_parameter("Model's scale_y cannot be negative or zero");
	}
	// When adjust_mask=false we check the mask's dimensionality later
	if (mask && adjust_mask && mask.getDimensions() != requested_dimensions) {
		throw invalid_parameter("Mask dimensions != model dimensions");
	}

	input_analysis analysis;
	analysis.convolution_required = std::any_of(profiles.begin(), profiles.end(),
	                                            std::mem_fn(&Profile::do_convolve));

	if (analysis.convolution_required && !psf) {
		throw invalid_parameter("No psf provided but profile(s) requested convolution");
	}

	/* Validate all profiles. Each profile can fail during validation */
	for(auto &profile: this->profiles) {
		profile->validate();
	}

	analyze_expansion_requirements(requested_dimensions, mask, psf,
	                               finesampling, analysis, adjust_mask);
	return analysis;
}

bool Model::needs_adjustment(const Mask &mask, unsigned int finesampling,
    const input_analysis &analysis)
{
	if (!mask) {
		return false;
	}
	return analysis.mask_needs_convolution || analysis.mask_needs_psf_padding || finesampling > 1;
}

void Model::adjust(Mask &mask, const Image &psf, unsigned int finesampling,
    const input_analysis &analysis)
{
	if (!mask) {
		return;
	}
	if (finesampling > 1) {
		mask = mask.upsample(finesampling);
	}
	if (analysis.mask_needs_psf_padding) {
		mask = mask.extend(analysis.drawing_dims, analysis.psf_padding);
	}
	if (analysis.mask_needs_convolution) {
		mask = mask.expand_by(psf.getDimensions() / 2);
	}
	assert(mask.getDimensions() == analysis.drawing_dims);
}

void Model::adjust(Mask &mask, const Dimensions &dims, const Image &psf,
    unsigned int finesampling)
{
	input_analysis analysis;
	analysis.convolution_required = psf.size() > 0;
	analyze_expansion_requirements(dims, mask, psf, finesampling, analysis, true);
	if (analysis.mask_needs_adjustment) {
		adjust(mask, psf, finesampling, analysis);
	}
}

Image Model::evaluate(Point &offset_out)
{
	auto analysis = analyze_inputs();

	/* so long folks! */
	if (dry_run) {
		inform_offset({0, 0}, offset_out);
		return Image{analysis.drawing_dims};
	}

	// Adjust mask before passing it down to profiles
	Point offset;
	Image image;
	Mask adjusted_mask;
	if (adjust_mask && analysis.mask_needs_adjustment) {
		adjusted_mask = mask;
		adjust(adjusted_mask, psf, finesampling, analysis);
		image = produce_image(adjusted_mask, analysis, offset);
	}
	else {
		if (!adjust_mask && mask.getDimensions() != analysis.drawing_dims) {
			std::ostringstream os;
			os << "Mask dimensions != drawing dimensions: "
			   << mask.getDimensions() << " != " << analysis.drawing_dims;
			throw invalid_parameter(os.str());
		}
		image = produce_image(mask, analysis, offset);
	}

	// Remove PSF padding if one was added, and downsample if necessary
	if (analysis.psf_padding) {
		auto crop_offset = analysis.psf_padding;
		auto crop_dims = analysis.drawing_dims - analysis.psf_padding * 2;
		if (!crop) {
			// We need to remove the padding effects from the uncropped
			// area. For that we see how much more extra padding was added
			// only due to the extra padding
			auto conv_actual_padding = convolver->padding(analysis.drawing_dims, psf.getDimensions());
			auto conv_intended_padding = convolver->padding(analysis.drawing_dims - analysis.psf_padding * 2, psf.getDimensions());
			auto offset_diff = conv_actual_padding.first - conv_intended_padding.first;
			auto dim_diff = conv_actual_padding.second - conv_intended_padding.second;
			crop_dims = image.getDimensions() - analysis.psf_padding * 2;
			crop_dims -= offset_diff + dim_diff;
			crop_offset += offset_diff;
			offset -= offset_diff;
		}
		image = image.crop(crop_dims, crop_offset);
	}

	if (finesampling > 1 && !return_finesampled) {
		image = image.downsample(finesampling, Image::DownsamplingMode::SUM);
		offset /= finesampling;
	}
	// Only in this case we know exactly what to mask out; otherwise
	// users should have the original mask
	if (adjusted_mask) {
		image &= mask;
	}

	inform_offset(offset, offset_out);
	return image;
}

Image Model::produce_image(const Mask &mask, const input_analysis &analysis,
    Point &offset)
{
	// Avoiding memory allocation if no convolution is needed
	Image model_image{analysis.drawing_dims};
	Image to_convolve;
	if (analysis.convolution_required) {
		to_convolve = Image{analysis.drawing_dims};
	}

	// Capturing as a lambda to destroy "scratch" quickly
	[this, &model_image, &to_convolve, &analysis, &mask]() {
		Image scratch{analysis.drawing_dims};
		for(auto &profile: this->profiles) {
			profile->adjust_for_finesampling(finesampling);
			std::fill(scratch.begin(), scratch.end(), 0);
			profile->evaluate(scratch, mask,
				{scale.first / finesampling, scale.second / finesampling},
				analysis.psf_padding, magzero);
			if (profile->do_convolve()) {
				to_convolve += scratch;
			}
			else {
				model_image += scratch;
			}
		}
	}();

	// Perform convolution if needed, then add back to the model image
	offset = {0, 0};
	if (analysis.convolution_required) {
		Image psf_img(psf);
		psf_img.normalize();
		to_convolve = ensure_convolver()->convolve(to_convolve, psf_img, mask, crop, offset);
		// The result of the convolution might be bigger that the original,
		// dependingo on user settings, so we need to account for that
		if (to_convolve.getDimensions() != analysis.drawing_dims) {
			model_image = model_image.extend(to_convolve.getDimensions(), offset);
		}
		model_image += to_convolve;
	}

	/* Done! Good job :-) */
	return model_image;
}

std::map<std::string, std::shared_ptr<ProfileStats>> Model::get_stats() const {
	std::map<std::string, std::shared_ptr<ProfileStats>> stats;
	for(auto &p: profiles) {
		stats[p->get_name()] = p->get_stats();
	}
	return stats;
}

#ifdef PROFIT_DEBUG
std::map<std::string, std::map<int, int>> Model::get_profile_integrations() const {
	std::map<std::string, std::map<int, int>> profile_integrations;
	for(auto &p: profiles) {
		RadialProfile *rp = dynamic_cast<RadialProfile *>(p.get());
		if( rp ) {
			profile_integrations[rp->get_name()] = rp->get_integrations();
		}
	}
	return profile_integrations;
}
#endif /* PROFIT_DEBUG */

} /* namespace profit */