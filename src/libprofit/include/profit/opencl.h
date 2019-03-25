/**
 * User-facing header file for OpenCL functionality
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
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

#ifndef PROFIT_OPENCL_H
#define PROFIT_OPENCL_H

#include <map>
#include <memory>
#include <string>

#include "profit/config.h"
#include "profit/common.h"

namespace profit
{

/**
 * A datatype for storing an OpenCL version.
 * It should have the form major*100 + minor*10 (e.g., 120 for OpenCL 1.2)
 */
typedef unsigned int cl_ver_t;

/**
 * A structure holding two times associated with OpenCL commands:
 * submission and execution
 */
struct PROFIT_API OpenCL_command_times {
	OpenCL_command_times();
	nsecs_t submit;
	nsecs_t exec;
	OpenCL_command_times &operator+=(const OpenCL_command_times &other);
	const OpenCL_command_times operator+(const OpenCL_command_times &other) const;
};

/**
 * A structure holding a number of OpenCL command times (filling, writing,
 * kernel and reading) plus other OpenCL-related times.
 */
struct PROFIT_API OpenCL_times {
	OpenCL_times();
	nsecs_t kernel_prep;
	unsigned int nwork_items;
	OpenCL_command_times writing_times;
	OpenCL_command_times reading_times;
	OpenCL_command_times filling_times;
	OpenCL_command_times kernel_times;
	nsecs_t total;
};

/**
 * An OpenCL environment.
 *
 * This class holds all the required information to make libprofit work
 * against a given device in a particular platform.
 */
class PROFIT_API OpenCLEnv {

public:
	virtual ~OpenCLEnv() {};

	/**
	 * Returns the maximum OpenCL version supported by the underlying device.
	 */
	virtual cl_ver_t get_version() = 0;

	/**
	 * Returns the name of the OpenCL platform of this environment.
	 *
	 * @return The name of the OpenCL platform
	 */
	virtual std::string get_platform_name() = 0;
	/**
	 * Returns the name of the OpenCL device of this environment.
	 * @return The name of the OpenCL device
	 */
	virtual std::string get_device_name() = 0;

};

/// Handy typedef for shared pointers to OpenCL_env objects
/// This is the object users handle and pass around.
typedef std::shared_ptr<OpenCLEnv> OpenCLEnvPtr;

/**
 * A structure holding information about a specific OpenCL device
 */
typedef struct PROFIT_API _OpenCL_dev_info {

	/** The name of the device */
	std::string name;

	/** The OpenCL version supported by this device */
	cl_ver_t cl_version;

	/** Whether or not this device supports double floating-point precision */
	bool double_support;

} OpenCL_dev_info;

/**
 * An structure holding information about a specific OpenCL platform.
 */
typedef struct PROFIT_API _OpenCL_plat_info {

	/** The name of the platform */
	std::string name;

	/** The supported OpenCL version */
	cl_ver_t supported_opencl_version;

	/** A map containing information about all devices on this platform */

	std::map<int, OpenCL_dev_info> dev_info;
} OpenCL_plat_info;

/**
 * Queries the system about the OpenCL supported platforms and devices and returns
 * the information the caller.
 *
 * @return A map keyed by index, containing the information of each of the
 *         OpenCL platforms found on this system.
 */
PROFIT_API std::map<int, OpenCL_plat_info> get_opencl_info();

/**
 * Prepares an OpenCL working space for using with libprofit.
 *
 * This method will get the requested device on the requested platform, compile
 * the libprofit OpenCL kernel sources to be used against it, and set up a queue
 * on the device.
 *
 * @param platform_idx The index of the platform to use
 * @param device_idx The index of device to use in the platform
 * @param use_double Whether double floating-point support should be used in
 *        the device or not.
 * @param enable_profiling Whether OpenCL profiling capabilities should be
 *        turned on in the OpenCL Queue created within this envinronment.
 * @return A pointer to a OpenCL_env structure, which contains the whole set of
 *         elements required to work with the requested device.
 */
PROFIT_API OpenCLEnvPtr get_opencl_environment(
	unsigned int platform_idx,
	unsigned int device_idx,
	bool use_double,
	bool enable_profiling);

} /* namespace profit */

#endif /* PROFIT_OPENCL_H */