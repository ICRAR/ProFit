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


using namespace std;

namespace profit {

invalid_parameter::invalid_parameter(const string &what_arg) :
	exception(),
	m_what(what_arg)
{
	// no-op
}

invalid_parameter::~invalid_parameter() throw () {
	// no-op
}

const char *invalid_parameter::what() const throw() {
	return m_what.c_str();
}

#ifdef PROFIT_OPENCL
opencl_error::opencl_error(const string &what_arg) :
	exception(),
	m_what(what_arg)
{
	// no-op
}

opencl_error::~opencl_error() throw () {
	// no-op
}

const char *opencl_error::what() const throw() {
	return m_what.c_str();
}
#endif /* PROFIT_OPENCL */

} /* namespace profit */