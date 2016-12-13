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
#include <typeinfo>

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/profile.h"


using namespace std;

namespace profit {

Profile::Profile(const Model &model, const string &name) :
	model(model),
	name(name),
	convolve(false)
{
	// no-op
}

Profile::~Profile()
{
	// no-op
}

bool Profile::do_convolve() const {
	return convolve;
}

const string& Profile::get_name() const {
	return name;
}

template <typename T>
void Profile::set_parameter(const string &name, T val) {
	if( !parameter_impl(name, val) ) {
		ostringstream os;
		os << "Unknown " << typeid(val).name() << " parameter '" << name << "'";
		throw invalid_parameter(os.str());
	}
}

void Profile::parameter(const string &name, bool val) {
	set_parameter<bool>(name, val);
}

void Profile::parameter(const string &name, double val) {
	set_parameter<double>(name, val);
}

void Profile::parameter(const string &name, unsigned int val) {
	set_parameter<unsigned int>(name, val);
}

bool Profile::parameter_impl(const string &name, bool val) {

	if( name == "convolve" ) {
		convolve = val;
		return true;
	}

	return false;
}

bool Profile::parameter_impl(const string &name, double val) {
	return false;
}

bool Profile::parameter_impl(const string &name, unsigned int val) {
	return false;
}

} /* namespace profit */