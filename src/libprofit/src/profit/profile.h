/**
 * Header file for the base Profile class of libprofit
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

#ifndef PROFIT_PROFILE_H
#define PROFIT_PROFILE_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "profit/config.h"
#include "profit/common.h"
#include "profit/image.h"
#include "profit/opencl.h"

namespace profit
{

/* Forward declaration */
class Model;
class PixelScale;

/// Statistics for profile evaluations
class PROFIT_API ProfileStats {
public:
	// We need at least one virtual function so we create a polymorphic hierarchy
	// rooted at this class. Otherwise we cannot distinguish at runtime between
	// this and inherited classes via dynamic_cast, which we do in profit-cli.
	// An alternative would be to add a tag or enumeration here and have
	// inheriting subclasses use different values (a la C)
	virtual ~ProfileStats() {};
	nsecs_t total;
};

/// Subsampling statistics for radial profile evaluations
struct PROFIT_API radial_subsampling_stats {
	nsecs_t pre_subsampling;
	nsecs_t new_subsampling;
	nsecs_t inital_transform;
	OpenCL_times cl_times;
	nsecs_t final_transform;
	nsecs_t total;
};

/// Statistics for radial profile evaluations
struct PROFIT_API RadialProfileStats : ProfileStats {
	OpenCL_times cl_times;
	radial_subsampling_stats subsampling;
	nsecs_t final_image;
};

/**
 * The base profile class
 */
class PROFIT_API Profile {

public:

#if __cpp_alias_templates == 200704
	template <typename T>
	using parameter_holder = std::map<std::string, std::reference_wrapper<T>>;
#else
	template <typename T>
	class parameter_holder : public std::map<std::string, std::reference_wrapper<T>> { };
#endif

	/**
	 * Constructor
	 *
	 * @param model The model this profile belongs to
	 * @param name The name of this profile
	 */
	Profile(const Model & model, const std::string &name);

	/**
	 * Destructor
	 */
	virtual ~Profile() = default;

	/**
	 * Performs the initial profile validation, making sure that all parameters
	 * of the profile are correct and can be safely used to create an image.
	 * This function can signal an error by throwing an invalid_parameter exception.
	 */
	virtual void validate() = 0;

	/**
	 * Adjusts the internal parameters of this profile to the given finesampling parameter.
	 * Finesampling produces bigger images, and therefore any parameters that indicate a position
	 * in the image coordinate space needs to be adjusted (by multiplying it by the finesampling
	 * factor).
	 *
	 * @param finesampling The finesampling factor.
	 */
	virtual void adjust_for_finesampling(unsigned int finesampling);

	/**
	 * Performs the profile evaluation and saves the resulting image into
	 * the given @p image. If @p mask is not empty it has the same
	 * dimensions than the image, and only the pixels where the mask is set
	 * need to be calculated. The image's pixel is scale is given by @p scale,
	 * and the zero magnitude value is given by @p magzero (not used by all profiles,
	 * but by most).
	 *
	 * This is the main function of the profile.
	 *
	 * @param image The Image object where values need to be stored.
	 * @param mask The mask to apply during profile calculation.
	 * @param scale The pixel scale of the image.
	 * @param offset The offset of the profile's origin with respect to the
	 * the image's origin
	 * @param magzero The profile's zero magnitude value.
	 */
	virtual void evaluate(Image &image, const Mask &mask, const PixelScale &scale,
	    const Point &offset, double magzero) = 0;

	/**
	 * Parses @p parameter_spec, which should look like `name = value`, and
	 * sets that parameter value on the profile.
	 * @param parameter_spec The parameter name
	 * @throws invalid_parameter if @p parameter_spec fails to parse, or the
	 * parameter's value cannot be parsed correctly
	 * @throws unknown_parameter if @p parameter_spec refers to a parameter not
	 * supported by this profile
	 */
	void parameter(const std::string &parameter_spec);

	/**
	 * Sets the parameter `name` to `value`.
	 * @param name The parameter name
	 * @param value The parameter value
	 * @throws invalid_parameter if `name` corresponds with no known parameter
	 * on this profile of type `bool`.
	 */
	void parameter(const std::string &name, bool value);

	/**
	 * Sets the parameter `name` to `value`.
	 * @param name The parameter name
	 * @param value The parameter value
	 * @throws invalid_parameter if `name` corresponds with no known parameter
	 * on this profile of type `double`.
	 */
	void parameter(const std::string &name, double value);

	/**
	 * Sets the parameter `name` to `value`.
	 * @param name The parameter name
	 * @param value The parameter value
	 * @throws invalid_parameter if `name` corresponds with no known parameter
	 * on this profile of type `unsigned int`.
	 */
	void parameter(const std::string &name, unsigned int value);

	/**
	 * Returns the name of this profile
	 *
	 * @return the name of this profile
	 */
	const std::string& get_name(void) const;

	/**
	 * Returns whether the image generated by this profile needs to be
	 * convolved or not.
	 *
	 * @return whether to convolve (or not) the result of evaluating this profile
	 */
	bool do_convolve(void) const;

	/**
	 * Returns the runtime statistics after evaluating this profile.
	 * @return
	 */
	std::shared_ptr<ProfileStats> get_stats() const;

protected:

	/**
	 * Registers a boolean variable as a profile parameter. It is meant to be
	 * called from the different profiles' constructors to register their
	 * private members holding boolean parameter values.
	 *
	 * @param name The name of the boolean parameter
	 * @param variable The boolean variable holding the parameter
	 */
	void register_parameter(const char *name, bool &variable);

	/**
	 * Like register_parameter(const char *, bool), but for `unsigned int`
	 * parameters
	 *
	 * @param name The name of the unsigned int parameter
	 * @param variable The unsigned int variable holding the parameter
	 * @see register_parameter(const char *, bool)
	 */
	void register_parameter(const char *name, unsigned int &variable);

	/**
	 * Like register_parameter(const char *, bool), but for `double`
	 * parameters
	 *
	 * @param name The name of the double parameter
	 * @param variable The double variable holding the parameter
	 * @see register_parameter(const char *, bool)
	 */
	void register_parameter(const char *name, double &variable);

	/**
	 * A (constant) reference to the model this profile belongs to
	 */
	const Model &model;

	/**
	 * The name of this profile
	 */
	const std::string name;

private:

	/** @name Profile Parameters */
	// @{
	/**
	 * Whether the resulting image of this profile should be convolved or not.
	 */
	bool convolve;
	// @}

	parameter_holder<bool> bool_parameters;
	parameter_holder<unsigned int> uint_parameters;
	parameter_holder<double> double_parameters;

	std::shared_ptr<ProfileStats> stats;

	// RadialProfile sets a different type of stats, and until we have a more
	// generic stats API we simply let it use our private member
	friend class RadialProfile;
};

/// A pointer to a Profile object
typedef std::shared_ptr<Profile> ProfilePtr;

} /* namespace profit */

#endif /* PROFIT_PROFILE_H */
