/**
 * Base Profile class implementation
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

#include <sstream>
#include <string>

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/profile.h"
#include "profit/utils.h"


namespace profit {

Profile::Profile(const Model &model, const std::string &name) :
	model(model),
	name(name),
	convolve(false),
	stats()
{
	register_parameter("convolve", convolve);
}

void Profile::adjust_for_finesampling(unsigned int finesampling)
{
	// no-op
}

bool Profile::do_convolve() const {
	return convolve;
}

const std::string& Profile::get_name() const {
	return name;
}

std::shared_ptr<ProfileStats> Profile::get_stats() const {
	return stats;
}

void Profile::register_parameter(const char *name, bool &parameter)
{
	bool_parameters.insert({name, parameter});
}

void Profile::register_parameter(const char *name, unsigned int &parameter)
{
	uint_parameters.insert({name, parameter});
}

void Profile::register_parameter(const char *name, double &parameter)
{
	double_parameters.insert({name, parameter});
}

template <typename T>
void set_parameter(
	Profile::parameter_holder<T> &parameters,
	const std::string &name,
	const std::string &profile_name,
	T val)
{
	constexpr auto tname = type_info<T>::name;
	if (parameters.find(name) == parameters.end()) {
		std::ostringstream os;
		os << "Unknown " << tname << " parameter in profile " << profile_name << ": " << name;
		throw invalid_parameter(os.str());
	}
	parameters.at(name).get() = val;
}

template <typename T, typename Converter>
bool set_parameter(
	Profile::parameter_holder<T> &parameters,
	const std::string &name,
	const std::string &profile_name,
	const std::string &val,
	Converter &&converter)
{
	if (parameters.find(name) == parameters.end()) {
		return false;
	}

	try {
		T bval = converter(val);
		parameters.at(name).get() = bval;
		return true;
	} catch (const std::invalid_argument &e) {
		UNUSED(e);
		constexpr auto type_name = type_info<T>::name;
		std::ostringstream os;
		os << "Parameter " << name << " in profile " << profile_name;
		os << " is of type " << type_name << ", but given value cannot be parsed";
		os << " as " << type_name << ": " << val;
		throw invalid_parameter(os.str());
	}
}

void Profile::parameter(const std::string &name, bool val) {
	set_parameter(bool_parameters, name, get_name(), val);
}

void Profile::parameter(const std::string &name, double val) {
	set_parameter(double_parameters, name, get_name(), val);
}

void Profile::parameter(const std::string &name, unsigned int val) {
	set_parameter(uint_parameters, name, get_name(), val);
}

void Profile::parameter(const std::string &param_spec)
{
	auto parts = split(param_spec, "=");
	if (parts.size() != 2) {
		std::ostringstream os;
		os << "missing = in parameter: " << param_spec;
		throw invalid_parameter(os.str());
	}

	auto &pname = trim(parts[0]);
	auto &val = trim(parts[1]);

	bool found = (
		set_parameter(bool_parameters, pname, get_name(), val, [](const std::string &s) { return std::stoul(s, nullptr, 10); }) ||
		set_parameter(uint_parameters, pname, get_name(), val, [](const std::string &s) { return stoui(s); }) ||
		set_parameter(double_parameters, pname, get_name(), val, [](const std::string &s) { return std::stod(s); })
	);

	if (!found) {
		std::ostringstream os;
		os << "Profile " << get_name() << " doesn't support parameter " << pname;
		os << ", or parameter has invalid value: " << val;
		throw unknown_parameter(os.str());
	}
}

} /* namespace profit */