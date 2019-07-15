/**
 * Exception classes implementation
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

#include <string>

#include "profit/common.h"
#include "profit/exceptions.h"


namespace profit {

exception::exception(const std::string &what_arg) :
	m_what(what_arg)
{
	// no-op
}

exception::~exception() throw () {
	// no-op
}

const char *exception::what() const throw() {
	return m_what.c_str();
}

invalid_parameter::invalid_parameter(const std::string &what_arg) :
	exception(what_arg)
{
	// no-op
}

invalid_parameter::~invalid_parameter() throw () {
	// no-op
}

unknown_parameter::unknown_parameter(const std::string &what_arg) :
	invalid_parameter(what_arg)
{
	// no-op
}

unknown_parameter::~unknown_parameter() throw () {
	// no-op
}

opencl_error::opencl_error(const std::string &what_arg) :
	exception(what_arg)
{
	// no-op
}

opencl_error::~opencl_error() throw () {
	// no-op
}

fft_error::fft_error(const std::string &what_arg) :
	exception(what_arg)
{
	// no-op
}

fft_error::~fft_error() throw () {
	// no-op
}

fs_error::fs_error(const std::string &what_arg) :
	exception(what_arg)
{
	// no-op
}

fs_error::~fs_error() throw () {
	// no-op
}
} /* namespace profit */