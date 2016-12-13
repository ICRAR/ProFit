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

#include <string>
#include <vector>

namespace profit
{

/* Forward declaration */
class Model;

/**
 * The base profile class
 */
class Profile {

public:

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
	virtual ~Profile() = 0;

	/**
	 * Performs the initial profile validation, making sure that all parameters
	 * of the profile are correct and can be safely used to create an image.
	 * This function can signal an error by throwing an invalid_parameter exception.
	 */
	virtual void validate() = 0;

	/**
	 * Performs the profile evaluation and saves the resulting image into
	 * the given `image` array. This is the main function of the profile.
	 *
	 * @param image The vector where image values need to be stored.
	 *              Its size is `model.width` * `model.height`.
	 *              The data is organized by rows first, columns later;
	 *               i.e pixel (x,y) is accessed by `image[y*width + x]`
	 */
	virtual void evaluate(std::vector<double> &image) = 0;

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

	const std::string& get_name(void) const;

	bool do_convolve(void) const;

protected:

	/**
	 * Sets the parameter `name` to `value`. This method is meant to be
	 * overwritten by classes, and therefore care should be taken to check
	 * the return value from the parent class's implementation.
	 *
	 * @param name The parameter name
	 * @param value The parameter value
	 * @return Whether the parameter was set or not.
	 */
	virtual bool parameter_impl(const std::string &name, bool value);

	/**
	 * @see parameter_impl(const std::string, bool)
	 */
	virtual bool parameter_impl(const std::string &name, double value);

	/**
	 * @see parameter_impl(const std::string, bool)
	 */
	virtual bool parameter_impl(const std::string &name, unsigned int value);

	/**
	 * A (constant) reference to the model this profile belongs to
	 */
	const Model &model;

	/**
	 * The name of this profile
	 */
	const std::string name;

	/** @name Profile Parameters */
	// @{
	/**
	 * Whether the resulting image of this profile should be convolved or not.
	 */
	bool convolve;
	// @}

private:

	template <typename T>
	void set_parameter(const std::string &name, T value);

};

} /* namespace profit */

#endif /* PROFIT_PROFILE_H */
